// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_POLY_HASH_TABLE_GUARD
#define MATHICGB_POLY_HASH_TABLE_GUARD

#include "PolyRing.hpp"
#include "Poly.hpp"
#include <utility>
#include <memtailor.h>
#include <vector>

MATHICGB_NAMESPACE_BEGIN

class PolyHashTable {
public:
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::ConstMonoRef ConstMonoRef;
  typedef Monoid::MonoRef MonoRef;
  typedef coefficient Value;

  class Node {
  public:
    ConstMonoRef mono() const {return *Monoid::toMonoPtr(mMono);}
    MonoRef mono() {return *Monoid::toMonoPtr(mMono);}

    Value& value() {return mValue;}
    const Value& value() const {return mValue;}

  private:
    friend class PolyHashTable;

    Node*& next() {return mNext;}
    Node* next() const {return mNext;}

    Node* mNext;
    Value mValue;
    exponent mMono[1];
  };

  // Construct a hash table with at least requestedBucketCount buckets. There
  // may be more buckets. Currently the number is rounded up to the next power
  // of two.
  PolyHashTable(const PolyRing& ring):
    mHashToIndexMask(computeHashMask(1000)),
    mBuckets
      (make_unique_array<Node*>(hashMaskToBucketCount(mHashToIndexMask))),
    mRing(ring),
    mNodes(sizeofNode(ring)
    ),
    mSize()
  {
    mMaxSize = static_cast<size_t>(bucketCount() * maxLoadFactor());
    std::fill_n(mBuckets.get(), bucketCount(), nullptr);
  }

  const PolyRing::Monoid& monoid() const {return mRing.monoid();}

  /// Return how many buckets the hash table has.
  size_t bucketCount() const {
    return hashMaskToBucketCount(mHashToIndexMask);
  }

  /// Return the number of elements (not the number of buckets).
  size_t size() const {return mSize;}

  MATHICGB_INLINE
  std::pair<Node*, bool> insertProduct(ConstMonoRef a, ConstMonoRef b) {
    auto newNode = new (mNodes.alloc()) Node();
    monoid().multiply(a, b, newNode->mono());
    const auto abHash = monoid().hash(newNode->mono());
    auto& bucket = mBuckets[hashToIndex(abHash)];

    for (auto node = bucket; node != nullptr; node = node->next()) {
      if (abHash != monoid().hash(node->mono()))
        continue;
      if (monoid().equal(newNode->mono(), node->mono())) {
        mNodes.free(newNode);
        return std::make_pair(node, false); // found a*b.
      }
    }

    mRing.coefficientSet(newNode->value(), 0);
    newNode->next() = bucket;
    bucket = newNode;
    ++mSize;
    if (mSize >= mMaxSize)
      rehash(bucketCount() * 2);
    return std::make_pair(newNode, true); // inserted mono
  }

  MATHICGB_INLINE
  std::pair<Node*, bool> insertProduct
    (ConstMonoRef a, ConstMonoRef b, Value add)
  {
    auto p = insertProduct(a, b);
    mRing.coefficientAddTo(p.first->value(), add);
    return p;
  }

  MATHICGB_INLINE
  std::pair<Node*, bool> insertProduct(const_term a, const_term b)
  {
    Value prod;
    mRing.coefficientMult(a.coeff, b.coeff, prod);
    return insertProduct(a.monom, b.monom, prod);
  }

  MATHICGB_INLINE
  void remove(Node* nodeToRemove) {
    MATHICGB_ASSERT(nodeToRemove != 0);
    MATHICGB_ASSERT(mNodes.fromPool(nodeToRemove));
    const auto index = hashToIndex(monoid().hash(nodeToRemove->mono()));
    auto nodePtr = &mBuckets[index];
    while (*nodePtr != nodeToRemove) {
      MATHICGB_ASSERT(*nodePtr != nullptr);
      nodePtr = &(*nodePtr)->next();
    }
    *nodePtr = nodeToRemove->next();
    mNodes.free(nodeToRemove);
    --mSize;
  }

  /// Removes all elements and optimizes internal resources. This is
  /// fast if there are no elements, so if you know that there are no
  /// elements and that many operations have happened since the last clear,
  /// then call clear for better cache performance. If there is even one
  /// element, then this takes linear time in the number of buckets.
  void clear() {
    if (!empty()) {
      std::fill_n(mBuckets.get(), bucketCount(), nullptr);
      mSize = 0;
    }
    mNodes.freeAllBuffers();
  }

  bool empty() const {return mSize == 0;}

  size_t getMemoryUse() const {
    return bucketCount() * sizeof(mBuckets[0]) + mNodes.getMemoryUse();
  }

private:
  void rehash(const size_t requestedBucketCount) {
    const auto newHashToIndexMask = computeHashMask(requestedBucketCount);
    const auto newBucketCount = hashMaskToBucketCount(newHashToIndexMask);
    auto newBuckets = make_unique_array<Node*>(newBucketCount);
    std::fill_n(newBuckets.get(), newBucketCount, nullptr);

    const auto bucketsEnd = mBuckets.get() + bucketCount();
    for (auto bucket = mBuckets.get(); bucket != bucketsEnd; ++bucket) {
      for (auto node = *bucket; node != 0;) {
        const auto hash = monoid().hash(node->mono());
        const auto newIndex = hashToIndex(hash, newHashToIndexMask);
        const auto next = node->next();
        node->next() = newBuckets[newIndex];
        newBuckets[newIndex] = node;
        node = next;
      }
    }

    mHashToIndexMask = newHashToIndexMask;
    mBuckets = std::move(newBuckets);
    mMaxSize = static_cast<size_t>(bucketCount() * maxLoadFactor());
  }

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
    return
      sizeof(Node) +
      sizeof(Value) -
      sizeof(exponent) +
      ring.maxMonomialByteSize();
  }

  /// The maximum allowed value of elementCount() / bucketCount().
  static double maxLoadFactor() {return 0.10;}

  size_t hashToIndex(const HashValue hash) const {
    const auto index = hashToIndex(hash, mHashToIndexMask);
    MATHICGB_ASSERT(index == hash % bucketCount());
    return index;
  }

  static size_t hashToIndex(const HashValue hash, const HashValue mask) {
    return hash & mask;
  }

  HashValue mHashToIndexMask;
  std::unique_ptr<Node*[]> mBuckets;
  const PolyRing& mRing;
  memt::BufferPool mNodes;
  size_t mSize;
  size_t mMaxSize;
};

/*
// The hash table is a map:  monomial => coeff
// Operations required on monomials:
//  hash (this will currently pick out one entry of a monomial)
//
// Operations required on coefficients:
//  add, maybe multiply too?
//  isZero?
//
// Does not take ownership of any of the monomials.
class PolyHashTable {
public:
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::MonoRef MonoRef;
  typedef Monoid::ConstMonoRef ConstMonoRef;
  typedef Monoid::MonoPtr MonoPtr;
  typedef Monoid::ConstMonoPtr ConstMonoPtr;
  typedef coefficient Value;

  class Node {
  public:
    const const_monomial& mono() {return mMonom;}
    const const_monomial& mono() const {return mMonom;}

    Value& value() {return coeff;}
    const Value& value() const {return coeff;}

  private:
    friend class PolyHashTable;

    Node*& next() {return mNext;}
    Node* next() const {return mNext;}

    Node* mNext;
    coefficient coeff;
    const_monomial mMonom;
  };

  PolyHashTable(const PolyRing& ring);


  const Monoid& monoid() const {return mRing.monoid();}

  std::pair<Node*, bool> insertProduct
    (ConstMonoRef a, ConstMonoRef b);

  std::pair<Node*, bool> insertProduct
    (ConstMonoRef a, ConstMonoRef b, Value add)
  {
    auto p = insertProduct(a, b);
    mRing.coefficientAddTo(p.first->value(), add);
    return p;
  }

  std::pair<Node*, bool> insertProduct(const_term a, const_term b)
  {
    Value prod;
    mRing.coefficientMult(a.coeff, b.coeff, prod);
    return insertProduct(a.monom, b.monom, prod);
  }

  // Removes the node from the hash table.
  void remove(Node* n);

  void clear();  // Clear the table, and memory areas.

  bool empty() const {return mNodeCount == 0;}

  size_t getMemoryUse() const;

private:
  typedef std::vector<Node*> MonomialArray;
  void resize(size_t new_nbits);  // Don't change the nodes, table, but do recreate hashtable_

  Node* makeNode(const_monomial monom);
  bool lookup_and_insert(const_monomial m, Node *&result);

  const PolyRing& mRing;
  std::vector<Node*> mHashTable;
  HashValue mHashMask; // this is the number, in binary:  00001111...1, where
                    // the number of 1's is mLogTableSize

  memt::Arena mArena; // space for monomials represented in this class.  Also nodes??

  size_t mTableSize;
  size_t mLogTableSize; // number of bits in the table: mTableSize should be 2^mLogTableSize

  size_t mNodeCount;  // current number of nodes in the hash table
  size_t mBinCount;
  size_t mMaxCountBeforeRebuild;

  size_t mMonomialSize;
};
*/
MATHICGB_NAMESPACE_END
#endif
