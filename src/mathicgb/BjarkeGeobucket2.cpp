// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "BjarkeGeobucket2.hpp"

#include "TypicalReducer.hpp"
#include "PolyHashTable.hpp"
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

class MonoMap {
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
    friend class MonoMap;

    Node*& next() {return mNext;}
    Node* next() const {return mNext;}

    Node* mNext;
    Value mValue;
    exponent mMono[1];
  };

  // Construct a hash table with at least requestedBucketCount buckets. There
  // may be more buckets. Currently the number is rounded up to the next power
  // of two.
  MonoMap(
    const size_t requestedBucketCount,
    const PolyRing& ring
  ):
    mHashToIndexMask(computeHashMask(requestedBucketCount)),
    mBuckets
      (make_unique_array<Node*>(hashMaskToBucketCount(mHashToIndexMask))),
    mRing(ring),
    mNodes(sizeofNode(ring)
    ),
    mSize()
  {
    std::fill_n(mBuckets.get(), bucketCount(), nullptr);
  }

  const PolyRing::Monoid& monoid() const {return mRing.monoid();}

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
  }

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
    return std::make_pair(newNode, true); // inserted mono
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

private:
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
};


class BjarkeGeobucket2 : public TypicalReducer {
public:
  BjarkeGeobucket2(const PolyRing& ring):
    mRing(ring),
    mMap(10000, ring),
    mQueue2(QueueConfiguration2(ring.monoid()))
  {}

  virtual std::string description() const {return "bjarke geo buckets";}

  void insertTail(const_term multiplier, const Poly *g1) {
    MATHICGB_ASSERT(g1 != 0);
    MATHICGB_ASSERT(g1->termsAreInDescendingOrder());

    if (g1->nTerms() <= 1)
      return;

    mNodesTmp.clear();
    auto it = g1->begin();
    const auto end = g1->end();
    for (++it; it != end; ++it) {
      auto p = mMap.insertProduct(it.getMonomial(), multiplier.monom);
      coefficient prod;
      mRing.coefficientMult(it.getCoefficient(), multiplier.coeff, prod);
      mRing.coefficientAddTo(p.first->value(), prod);
      if (p.second)
        mNodesTmp.emplace_back(p.first);
    }
    if (!mNodesTmp.empty())
      mQueue2.push(mNodesTmp.begin(), mNodesTmp.end());
  }

  void insert(monomial multiplier, const Poly *g1) {
    MATHICGB_ASSERT(g1 != 0);
    MATHICGB_ASSERT(g1->termsAreInDescendingOrder());

    mNodesTmp.clear();
    const auto end = g1->end();
    for (auto it = g1->begin(); it != end; ++it) {
      auto p = mMap.insertProduct(it.getMonomial(), multiplier);
      mRing.coefficientAddTo(p.first->value(), it.getCoefficient());
      if (p.second)
        mNodesTmp.emplace_back(p.first);
    }
    if (!mNodesTmp.empty())
      mQueue2.push(mNodesTmp.begin(), mNodesTmp.end());
  }

  virtual bool leadTerm(const_term& result) {
    while (!mQueue2.empty()) {
      const auto node = mQueue2.top();
      if (node->value() != 0) {
        result.coeff = node->value();
        result.monom = Monoid::toOld(node->mono());
        return true;
      }
      mQueue2.pop();
      mMap.remove(node);
    }
    return false;
  }

  virtual void removeLeadTerm() {
    MATHICGB_ASSERT(!mQueue2.empty());
    const auto node = mQueue2.top();
    mQueue2.pop();
    mMap.remove(node);
  }

  virtual size_t getMemoryUse() const {
    size_t result = TypicalReducer::getMemoryUse();
    //result += mMap.getMemoryUse();
    return result;
  }


protected:
  void resetReducer() {
    const_term t;
    while (!mQueue2.empty()) {
      const auto node = mQueue2.top();
      mQueue2.pop();
      mMap.remove(node);
    }
    MATHICGB_ASSERT(mMap.empty());
    MATHICGB_ASSERT(mQueue2.empty());
    mMap.clear();
  }

private:
  class QueueConfiguration {
  public:
    typedef PolyRing::Monoid Monoid;

    QueueConfiguration(const Monoid& monoid):
      mMonoid(monoid), geoBase(4), minBucketSize(1) {}

    typedef PolyHashTable::node* Entry;

    typedef bool CompareResult;
    CompareResult compare(const Entry& a, const Entry& b) const {
      return mMonoid.lessThan(a->monom, b->monom);
    }
    bool cmpLessThan(CompareResult r) const {return r;}

    static const bool supportDeduplication = false;
    bool cmpEqual(CompareResult r) const {
      MATHICGB_ASSERT(false); // Not supposed to be used.
      return false;
    }
    Entry deduplicate(const Entry& a, const Entry& /* b */) const {
      MATHICGB_ASSERT(false); // Not supposed to be used.
      return a;
    }

    static const bool minBucketBinarySearch = true;
    static const bool trackFront = true;
    static const bool premerge = false;
    static const bool collectMax = false;
    static const mic::GeobucketBucketStorage bucketStorage =
      static_cast<mic::GeobucketBucketStorage>(1);
    static const size_t insertFactor = 1;

    const size_t geoBase;
    const size_t minBucketSize;

  private:
    const Monoid& mMonoid;
  };

  class QueueConfiguration2 {
  public:
    typedef PolyRing::Monoid Monoid;

    QueueConfiguration2(const Monoid& monoid):
      mMonoid(monoid), geoBase(4), minBucketSize(1) {}

    typedef MonoMap::Node* Entry;

    typedef bool CompareResult;
    CompareResult compare(const Entry& a, const Entry& b) const {
      return mMonoid.lessThan(a->mono(), b->mono());
    }
    bool cmpLessThan(CompareResult r) const {return r;}

    static const bool supportDeduplication = false;
    bool cmpEqual(CompareResult r) const {
      MATHICGB_ASSERT(false); // Not supposed to be used.
      return false;
    }
    Entry deduplicate(const Entry& a, const Entry& /* b */) const {
      MATHICGB_ASSERT(false); // Not supposed to be used.
      return a;
    }

    static const bool minBucketBinarySearch = true;
    static const bool trackFront = true;
    static const bool premerge = false;
    static const bool collectMax = false;
    static const mic::GeobucketBucketStorage bucketStorage =
      static_cast<mic::GeobucketBucketStorage>(1);
    static const size_t insertFactor = 1;

    const size_t geoBase;
    const size_t minBucketSize;

  private:
    const Monoid& mMonoid;
  };

  mutable std::vector<MonoMap::Node*> mNodesTmp;
  const PolyRing& mRing;
  MonoMap mMap;
  mic::Geobucket<QueueConfiguration2> mQueue2;
};

std::unique_ptr<TypicalReducer> makeBjarkeGeobucket2(const PolyRing& ring) {
  return make_unique<BjarkeGeobucket2>(ring);
}

MATHICGB_NAMESPACE_END
