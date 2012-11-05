#ifndef MATHICGB_FIXED_SIZE_MONOMIAL_MAP_GUARD
#define MATHICGB_FIXED_SIZE_MONOMIAL_MAP_GUARD

#include "Atomic.hpp"
#include "PolyRing.hpp"
#include <memtailor.h>
#include <limits>
#include <vector>
#include <algorithm>
#include <mutex>

/// Concurrent hashtable mapping from monomials to T with a fixed number of
/// buckets. Lookups are lockless while insertions grab a lock.
///
/// There is no limitation on the number of entries that can be inserted,
/// but performance will suffer if the ratio of elements to buckets gets
/// high.
///
/// You can insert new values but you cannot change the value that an
/// already-inserted value maps to. It is possible to clear the table
/// but this operation is not safe for concurrency.
template<class T>
class FixedSizeMonomialMap {
public:
  typedef T mapped_type;
  typedef std::pair<const_monomial, mapped_type> value_type;

  // Construct a hash table with at least requestedBucketCount buckets. There
  // may be more buckets. Currently the number is rounded up to the next power
  // of two.
  FixedSizeMonomialMap(
    const size_t requestedBucketCount,
    const PolyRing& ring
  ):
    mHashToIndexMask(computeHashMask(requestedBucketCount)),
    mBuckets(
      make_unique_array<Atomic<Node*>>(hashMaskToBucketCount(mHashToIndexMask))
    ),
    mRing(ring),
    mNodeAlloc(sizeofNode(ring))
  {
    // Calling new int[x] does not zero the array. std::atomic has a trivial
    // constructor so the same thing is true of new atomic[x]. Calling
    // new int[x]() is supposed to zero initialize but this apparently
    // does not work on GCC. So we have to fill the table with nulls
    // manually. This was wonderful to debug btw.
    // We can store relaxed as the constructor does not run concurrently.
    setTableEntriesToNullRelaxed();
  }

  /// Construct a hash table with at least requestedBucketCount buckets and
  /// insert the elements from the parameter map.
  ///
  /// The parameter map remains a valid object that can satisfy queries.
  /// However, it is an error to call non-const methods on the map other
  /// than the destructor. Also, insertions into *this may or may not
  /// be reflected for queries into map and some of the entries currently
  /// in map will disappear.
  FixedSizeMonomialMap(
    const size_t requestedBucketCount,
    FixedSizeMonomialMap<T>&& map
  ):
    mHashToIndexMask(computeHashMask(requestedBucketCount)),
    mBuckets(
      make_unique_array<Atomic<Node*>>(hashMaskToBucketCount(mHashToIndexMask))
    ),
    mRing(map.ring()),
    mNodeAlloc(std::move(map.mNodeAlloc))
  {
    // We can store relaxed as the constructor does not run concurrently.
    setTableEntriesToNullRelaxed();
    const auto tableEnd = map.mBuckets.get() + map.bucketCount();
    for (auto tableIt = map.mBuckets.get(); tableIt != tableEnd; ++tableIt) {
      for (Node* node = tableIt->load(); node != 0;) {
        const size_t index = hashToIndex(mRing.monomialHashValue(node->mono));
        Node* const next = node->next.load();
        node->next.store(mBuckets[index].load());
        mBuckets[index].store(node);
        node = next;
      }
    }
  }

  /// Return how many buckets the hash table has.
  size_t bucketCount() const {
    return hashMaskToBucketCount(mHashToIndexMask);
  }

  const PolyRing& ring() const {return mRing;}

  // Returns the value associated to mono or null if there is no such value.
  const mapped_type* find(const const_monomial mono) const {
    const HashValue monoHash = mRing.monomialHashValue(mono);
    const Node* node = bucketAtIndex(hashToIndex(monoHash));
    for (; node != 0; node = node->next.load(std::memory_order_consume)) {
      // To my surprise, it seems to be faster to comment out this branch.
      // I guess the hash table has too few collisions to make it worth it.
      //if (monoHash != mRing.monomialHashValue(node->mono))
      //  continue;
      if (mRing.monomialEqualHintTrue(mono, node->mono))
        return &node->value;
    }
    return 0;
  }

  // As find on the product a*b.
  MATHICGB_INLINE const mapped_type* findProduct(
    const const_monomial a,
    const const_monomial b
  ) const {
    const HashValue abHash = mRing.monomialHashOfProduct(a, b);
    const Node* node = bucketAtIndex(hashToIndex(abHash));
    for (; node != 0; node = node->next.load(std::memory_order_consume)) {
      // To my surprise, it seems to be faster to comment out this branch.
      // I guess the hash table has too few collisions to make it worth it.
      //if (abHash != mRing.monomialHashValue(node->mono))
      //  continue;
      if (mRing.monomialIsProductOfHintTrue(a, b, node->mono))
        return &node->value;
    }
    return 0;
  }

  /// As findProduct but looks for a1*b and a2*b at one time.
  MATHICGB_INLINE
  std::pair<const mapped_type*, const mapped_type*> findTwoProducts(
    const const_monomial a1,
    const const_monomial a2,
    const const_monomial b
  ) const {
    const HashValue a1bHash = mRing.monomialHashOfProduct(a1, b);
    const HashValue a2bHash = mRing.monomialHashOfProduct(a2, b);
    const Node* const node1 = bucketAtIndex(hashToIndex(a1bHash));
    const Node* const node2 = bucketAtIndex(hashToIndex(a2bHash));

    if (node1 != 0 && node2 != 0 && mRing.monomialIsTwoProductsOfHintTrue
      (a1, a2, b, node1->mono, node2->mono)
    )
      return std::make_pair(&node1->value, &node2->value);
    else
      return std::make_pair(findProduct(a1, b), findProduct(a2, b));
  }

  /// Makes value.first map to value.second unless value.first is already
  /// present in the map - in that case nothing is done. If p is the returned
  /// pair then *p.first is the value that value.first maps to after the insert
  /// and p.second is true if an insertion was performed. *p.first will not
  /// equal value.second if an insertion was not performed - unless the
  /// inserted value equals the already present value.
  std::pair<const mapped_type*, bool> insert(const value_type& value) {
    const std::lock_guard<std::mutex> lockGuard(mInsertionMutex);
    // find() loads buckets with memory_order_consume, so it may seem like
    // we need some extra synchronization to make sure that we have the
    // most up to date view of the bucket that value.first goes in -
    // what if a pending insert is hiding in a cache of some other processor
    // somewhere? We in fact have an up to date view of every bucket in the
    // the table already because inserts only happen while holding the
    // insertion mutex and by locking that mutex we have synchronized with
    // all threads that previously did insertions.
    {
      const mapped_type* const found = find(value.first);
      if (found != 0)
        return std::make_pair(found, false); // key already present
    }

    const auto node = static_cast<Node*>(mNodeAlloc.alloc());
    const size_t index = hashToIndex(mRing.monomialHashValue(value.first));
    // the constructor initializes the first field of node->mono, so
    // it has to be called before copying the monomial.
    new (node) Node(bucketAtIndex(index), value.second);
    {
      Monomial nodeTmp(node->mono);
      ring().monomialCopy(value.first, nodeTmp);
    }
    // we cannot store with memory_order_relaxed here because unlocking the
    // lock only synchronizes with threads who later grab the lock - it does
    // not synchronize with reading threads since they do not grab the lock.
    mBuckets[index].store(node, std::memory_order_release);
    return std::make_pair(&node->value, true); // successful insertion
  }

  /// This operation removes all entries from the table. This operation
  /// requires synchronization with and mutual exclusion from all other
  /// clients of *this - you need to supply this synchronization manually.
  void clearNonConcurrent() {
#ifdef MATHICGB_DEBUG
    // requires mutual exclusion from both readers and writers, but we can
    // only assert on this for the writers.
    if (mInsertionMutex.try_lock())
      mInsertionMutex.unlock();
    else {
      MATHICGB_ASSERT(false);
    }
#endif
    // we can store relaxed as the client supplies synchronization.
    setTableEntriesToNullRelaxed();

    // This is the reason that we cannot support this operation concurrently -
    // we have no way to know when it is safe to deallocate the monomials
    // since readers do no synchronization.
    mNodeAlloc.freeAllBuffers();
  }

private:
  void setTableEntriesToNullRelaxed() {
    const auto tableEnd = mBuckets.get() + bucketCount();
    for (auto tableIt = mBuckets.get(); tableIt != tableEnd; ++tableIt)
      tableIt->store(0, std::memory_order_relaxed);
  }

  struct Node {
    Node(Node* next, const mapped_type value): next(next), value(value) {}

    Atomic<Node*> next;
    const mapped_type value;
    exponent mono[1];
  };

  static HashValue computeHashMask(const size_t requestedBucketCount) {
    // round request up to nearest power of 2.
    size_t pow2 = 1;
    while (pow2 < requestedBucketCount && 2 * pow2 != 0)
      pow2 *= 2;
    MATHICGB_ASSERT(pow2 > 0 && (pow2 & (pow2 - 1)) == 0); // power of two

    // If casting to a hash value overflows, then we get the maximum
    // possible number of buckets based on the range of the hash
    // value type. Only unsigned overflow is defined, so we need
    // to assert that the hash type is unsigned.
    static_assert(!std::numeric_limits<HashValue>::is_signed, "");
    const auto hashToIndexMask = static_cast<HashValue>(pow2 - 1);
    MATHICGB_ASSERT(pow2 == hashMaskToBucketCount(hashToIndexMask));
    return hashToIndexMask;
  }

  static size_t hashMaskToBucketCount(const HashValue mask) {
    const auto count = static_cast<size_t>(mask) + 1u; // should be power of 2
    MATHICGB_ASSERT(count > 0 && (count & (count - 1)) == 0); 
    return count;
  }

  static size_t sizeofNode(const PolyRing& ring) {
    return sizeof(Node) - sizeof(exponent) + ring.maxMonomialByteSize();
  }

  size_t hashToIndex(HashValue hash) const {
    const auto index = hash & mHashToIndexMask;
    MATHICGB_ASSERT(index == hash % bucketCount());
    return index;
  }

  Node* bucketAtIndex(size_t index) {
    MATHICGB_ASSERT(index < bucketCount());
    return mBuckets[index].load(std::memory_order_consume);
  }

  const Node* bucketAtIndex(size_t index) const {
    MATHICGB_ASSERT(index < bucketCount());
    return mBuckets[index].load(std::memory_order_consume);
  }

  const HashValue mHashToIndexMask;
  std::unique_ptr<Atomic<Node*>[]> const mBuckets;
  const PolyRing& mRing;
  memt::BufferPool mNodeAlloc; // nodes are allocated from here.
  std::mutex mInsertionMutex;
};

#endif
