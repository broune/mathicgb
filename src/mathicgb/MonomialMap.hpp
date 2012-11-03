#ifndef MATHICGB_MONOMIAL_MAP_GUARD
#define MATHICGB_MONOMIAL_MAP_GUARD

#include "PolyRing.hpp"
#include <memtailor.h>
#include <limits>
#include <vector>
#include <algorithm>

/// A mapping from monomials to MapToType. The interface is the same as
/// for STL maps. The interface is not complete. Add methods you need when
/// you need it.
template<class MapToType>
class MonomialMap;

template<class MTT>
class FixedSizeMonomialMap {
public:
  typedef MTT mapped_type;
  typedef std::pair<const_monomial, mapped_type> value_type;

private:
  struct Node {
    Node* next;
    mapped_type value;
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

public:
  size_t bucketCount() const {
    return hashMaskToBucketCount(mHashToIndexMask);
  }

  FixedSizeMonomialMap(
    const size_t requestedBucketCount,
    const PolyRing& ring
  ):
    mHashToIndexMask(computeHashMask(requestedBucketCount)),
    mBuckets(new Node*[hashMaskToBucketCount(mHashToIndexMask)]),
    mRing(ring),
    mNodeAlloc(sizeofNode(ring)),
    mEntryCount(0)
  {
    std::fill_n(mBuckets, bucketCount(), static_cast<Node*>(0));
  }

  ~FixedSizeMonomialMap() {
    delete[] mBuckets;
  }

  /// map must not be mutated except to destruct it after this operation.
  /// Queries into map will continue to give either a correct result or
  /// a spurious miss. There will not be spurious hits.
  FixedSizeMonomialMap(
    const size_t requestedBucketCount,
    FixedSizeMonomialMap<MTT>&& map
  ):
    mHashToIndexMask(computeHashMask(requestedBucketCount)),
    mBuckets(new Node*[hashMaskToBucketCount(mHashToIndexMask)]),
    mRing(map.ring()),
    mNodeAlloc(std::move(map.mNodeAlloc)),
    mEntryCount(0)
  {
    std::fill_n(mBuckets, bucketCount(), static_cast<Node*>(0));
    const auto tableEnd = map.mBuckets + map.bucketCount();
    for (auto tableIt = map.mBuckets; tableIt != tableEnd; ++tableIt) {
      for (Node* node = *tableIt; node != 0;) {
        const size_t index = hashToIndex(mRing.monomialHashValue(node->mono));
        Node* const next = node->next;
        node->next = mBuckets[index];
        mBuckets[index] = node;
        node = next;
      }
    }
  }

  const PolyRing& ring() const {return mRing;}
  size_t size() const {return mEntryCount;}

  const mapped_type* find(const const_monomial mono) const {
    const HashValue monoHash = mRing.monomialHashValue(mono);
    const Node* node = bucketAtIndex(hashToIndex(monoHash));
    for (; node != 0; node = node->next) {
      // To my surprise, it seems to be faster to comment out this branch.
      // I guess the hash table has too few collisions to make it worth it.
      //if (monoHash != mRing.monomialHashValue(node->mono))
      //  continue;
      if (mRing.monomialEqualHintTrue(mono, node->mono))
        return &node->value;
    }
    return 0;
  }

  MATHICGB_INLINE const mapped_type* findProduct(
    const const_monomial a,
    const const_monomial b
  ) const {
    const HashValue abHash = mRing.monomialHashOfProduct(a, b);
    const Node* node = bucketAtIndex(hashToIndex(abHash));
    for (; node != 0; node = node->next) {
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

  void insert(const value_type& value) {
    Node* node = static_cast<Node*>(mNodeAlloc.alloc());
    const size_t index = hashToIndex(mRing.monomialHashValue(value.first));
    {
      Monomial nodeTmp(node->mono);
      ring().monomialCopy(value.first, nodeTmp);
    }
    new (&node->value) mapped_type(value.second);
    node->next = bucketAtIndex(index);
    mBuckets[index] = node;
    ++mEntryCount;
  }

  void clear() {
    std::fill_n(mBuckets, bucketCount(), static_cast<Node*>(0));
    mNodeAlloc.freeAllBuffers();
    mEntryCount = 0;
  }

private:
  size_t hashToIndex(HashValue hash) const {
    const auto index = hash & mHashToIndexMask;
    MATHICGB_ASSERT(index == hash % bucketCount());
    return index;
  }

  Node* bucketAtIndex(size_t index) {
    MATHICGB_ASSERT(index < bucketCount());
    return mBuckets[index];
  }

  const Node* bucketAtIndex(size_t index) const {
    MATHICGB_ASSERT(index < bucketCount());
    return mBuckets[index];
  }

  const HashValue mHashToIndexMask;
  Node** const mBuckets;
  const PolyRing& mRing;
  memt::BufferPool mNodeAlloc; // nodes are allocated from here.
  size_t mEntryCount;
};



namespace MonomialMapInternal {
  // The map type MapClass is defined here. This is here in
  // this namespace to avoid cluttering the class definition with
  // a lot of ifdef's.

  template<class MTT>
  class MapClass {
  public:
    typedef MapClass Map;
    typedef std::pair<const_monomial, MTT> value_type;
    typedef MTT mapped_type;

  private:
    struct Node {
      Node* next;
      mapped_type value;
      exponent mono[1];
    };

  public:
    MapClass(const PolyRing& ring, size_t bucketCount = 400000):
      mRing(ring),
      mTable(),
      mNodeAlloc(sizeof(Node) - sizeof(exponent) + ring.maxMonomialByteSize())
    {
      mGrowWhenThisManyEntries = 0;
      mTableSize = mTable.size();
      mEntryCount = 0;
      growTable();
    }

    const PolyRing& ring() const {return mRing;}

    MATHICGB_INLINE const mapped_type* findProduct(
      const const_monomial a,
      const const_monomial b
    ) const {
      return reader().findProduct(a, b);
    }

    /// As findProduct but looks for a1*b and a2*b at one time.
    MATHICGB_INLINE
    std::pair<const mapped_type*, const mapped_type*> findTwoProducts(
      const const_monomial a1,
      const const_monomial a2,
      const const_monomial b
    ) const {
      return reader().findTwoProducts(a1, a2, b);
    }

    void insert(const value_type& value) {
      Node* node = static_cast<Node*>(mNodeAlloc.alloc());
      const size_t index = hashToIndex(mRing.monomialHashValue(value.first));
      {
        Monomial nodeTmp(node->mono);
        ring().monomialCopy(value.first, nodeTmp);
      }
      new (&node->value) mapped_type(value.second);
      node->next = entry(index);
      mTable[index] = node;
      ++mEntryCount;
      MATHICGB_ASSERT(mEntryCount <= mGrowWhenThisManyEntries);
      if (mEntryCount == mGrowWhenThisManyEntries)
        growTable();
    }

    void clear() {
      std::fill(mTable.begin(), mTable.end(), static_cast<Node*>(0));
      mNodeAlloc.freeAllBuffers();
      mEntryCount = 0;
    }

    size_t size() const {return mEntryCount;}

    class MapReader {
    public:
      const mapped_type* find(const const_monomial mono) const {
        const HashValue monoHash = mRing.monomialHashValue(mono);
        const Node* node = entry(hashToIndex(monoHash));
        for (; node != 0; node = node->next) {
          // To my surprise, it seems to be faster to comment out this branch.
          // I guess the hash table has too few collisions to make it worth it.
          //if (monoHash != mRing.monomialHashValue(node->mono))
          //  continue;
          if (mRing.monomialEqualHintTrue(mono, node->mono))
            return &node->value;
        }
        return 0;
      }

      const PolyRing& ring() {return mRing;}

      size_t bucketCount() const {
        MATHICGB_ASSERT(mHashToIndexMask <
          std::numeric_limits<decltype(mHashToIndexMask)>::max());
        return mHashToIndexMask + 1;
      }

      MATHICGB_INLINE const mapped_type* findProduct(
        const const_monomial a,
        const const_monomial b
      ) const {
        const HashValue abHash = mRing.monomialHashOfProduct(a, b);
        const Node* node = entry(hashToIndex(abHash));
        for (; node != 0; node = node->next) {
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
        const Node* const node1 = entry(hashToIndex(a1bHash));
        const Node* const node2 = entry(hashToIndex(a2bHash));

        if (node1 != 0 && node2 != 0 && mRing.monomialIsTwoProductsOfHintTrue
          (a1, a2, b, node1->mono, node2->mono)
        )
          return std::make_pair(&node1->value, &node2->value);
        else
          return std::make_pair(findProduct(a1, b), findProduct(a2, b));
      }

    private:
      friend class MapClass<MTT>;
      size_t hashToIndex(const HashValue hash) const {
        const auto index = hash & mHashToIndexMask;
        MATHICGB_ASSERT(index == hash % bucketCount());
        return index;
      }

      const Node* entry(const size_t index) const {
        MATHICGB_ASSERT(index < bucketCount());
        return mTable[index];
      }

      MapReader(
        const PolyRing& ring,
        const Node* const* const table,
        const HashValue mHashToIndexMask
      ):
        mRing(ring),
        mTable(table),
        mHashToIndexMask(mHashToIndexMask) {}

      const PolyRing& mRing;
      const Node* const* const mTable;
      const HashValue mHashToIndexMask;
    };

    const mapped_type* find(const const_monomial mono) const {
      return reader().find(mono);
    }

    const MapReader reader() const {
      return MapReader(ring(), mTable.data(), mHashToIndexMask);
    }

  private:
    size_t hashToIndex(HashValue hash) const {
      const auto index = hash & mHashToIndexMask;
      MATHICGB_ASSERT(index == hash % mTable.size());
      return index;
    }

    Node* entry(size_t index) {
      MATHICGB_ASSERT(index < mTable.size());
      return mTable[index];
    }

    const Node* entry(size_t index) const {
      MATHICGB_ASSERT(index < mTable.size());
      return mTable[index];
    }

    void growTable() {
      // Determine parameters for larger hash table
      const size_t initialTableSize = 1 << 16; // must be a power of two!!!
      const float maxLoadFactor = 0.33f; // max value of size() / mTable.size()
      const float growthFactor = 2.0f; // multiply table size by this on growth

      const size_t newTableSize = mTable.empty() ?
        initialTableSize :
        static_cast<size_t>(mTable.size() * growthFactor + 0.5f); // round up
      const auto newGrowWhenThisManyEntries =
        static_cast<size_t>(newTableSize / maxLoadFactor); // round down

      MATHICGB_ASSERT((newTableSize & (newTableSize - 1)) == 0); // power of two

      // Move nodes from current table into new table
      decltype(mTable) newTable(newTableSize);
      HashValue newHashToIndexMask = static_cast<HashValue>(newTableSize - 1);
      const auto tableEnd = mTable.end();
      for (auto tableIt = mTable.begin(); tableIt != tableEnd; ++tableIt) {
        Node* node = *tableIt;
        while (node != 0) {
          const size_t index =
            mRing.monomialHashValue(node->mono) & newHashToIndexMask;
          MATHICGB_ASSERT(index < newTable.size());
          Node* const next = node->next;
          node->next = newTable[index];
          newTable[index] = node;
          node = next;
        }
      }

      // Move newly calculated table into place
      mTableSize = newTableSize;      
      mTable = std::move(newTable);
      mHashToIndexMask = newHashToIndexMask;
      mGrowWhenThisManyEntries = newGrowWhenThisManyEntries;

      MATHICGB_ASSERT(mTableSize < mGrowWhenThisManyEntries);
    }

    size_t mGrowWhenThisManyEntries;
    size_t mTableSize;
    HashValue mHashToIndexMask;
    size_t mEntryCount;
    const PolyRing& mRing;
    std::vector<Node*> mTable;
    memt::BufferPool mNodeAlloc; // nodes are allocated from here.
  };
}


template<class MTT>
class MonomialMap {
public:
  typedef MTT mapped_type;
  //typedef MonomialMapInternal::MapClass<MTT> MapType;
  typedef FixedSizeMonomialMap<MTT> MapType;

  //typedef typename MapType::Map::iterator iterator;
  //typedef typename MapType::Map::const_iterator const_iterator;
  typedef typename MapType::value_type value_type;

  /// Can be used for non-mutating accesses to the monomial map. This reader
  /// can miss queries that should hit if the map is mutated after the
  /// reader was constructed - it may only give access to a subset of the
  /// entries in this case. The only guarantee is that if the reader
  /// reports a hit then it really is a hit and spurious misses only happen
  /// after the map has been mutated.
  class SubsetReader {
  public:
    SubsetReader(const MonomialMap<MTT>& map): mMap(map.mMap.get()) {}

    const mapped_type* find(const_monomial m) const {
      return mMap->find(m);
    }

    const mapped_type* findProduct(
      const const_monomial a,
      const const_monomial b
    ) const {
      return mMap->findProduct(a, b);
    }

    MATHICGB_INLINE
    std::pair<const mapped_type*, const mapped_type*> findTwoProducts(
      const const_monomial a1,
      const const_monomial a2,
      const const_monomial b
    ) const {
      return mMap->findTwoProducts(a1, a2, b);
    }

  private:
    const FixedSizeMonomialMap<MTT>* const mMap;
  };

  MonomialMap(const PolyRing& ring):
    mMap(make_unique<MapType>(InitialBucketCount, ring)),
    mCapacityUntilGrowth(maxEntries(mMap->bucketCount())) {}

  const mapped_type* find(const_monomial m) const {
    return mMap->find(m);
  }

  const mapped_type* findProduct(
    const const_monomial a,
    const const_monomial b
  ) const {
    return mMap->findProduct(a, b);
  }

  MATHICGB_INLINE
  std::pair<const mapped_type*, const mapped_type*> findTwoProducts(
    const const_monomial a1,
    const const_monomial a2,
    const const_monomial b
  ) const {
    return mMap->findTwoProducts(a1, a2, b);
  }

  //size_t size() const {return mMap->size();}
  void clear() {mMap->clear();}

  void insert(const value_type& val) {
    while (mCapacityUntilGrowth == 0) {
      const size_t currentSize = size();
      if (mMap->bucketCount() >
        std::numeric_limits<size_t>::max() / GrowthFactor)
      {
        throw std::bad_alloc();
      }
      const size_t newBucketCount = mMap->bucketCount() * GrowthFactor;
      auto nextMap = make_unique<MapType>(newBucketCount, std::move(*mMap));
      mOldMaps.push_back(std::move(mMap));
      mMap = std::move(nextMap);
      mCapacityUntilGrowth = maxEntries(mMap->bucketCount()) - currentSize;
    }
    mMap->insert(val);
    MATHICGB_ASSERT(mCapacityUntilGrowth > 0);
    --mCapacityUntilGrowth;
  }

  size_t size() const {
    return maxEntries(mMap->bucketCount()) - mCapacityUntilGrowth;
  }

  const PolyRing& ring() const {return mMap->ring();}

private:
  static const size_t MinBucketsPerEntry = 3; // inverse of max load factor
  static const size_t GrowthFactor = 2;
  static const size_t InitialBucketCount = 1 << 16;

  static size_t maxEntries(const size_t bucketCount) {
    return (bucketCount + (MinBucketsPerEntry - 1)) / MinBucketsPerEntry;
  }

  std::unique_ptr<MapType> mMap;
  size_t mCapacityUntilGrowth;
  std::vector<std::unique_ptr<MapType>> mOldMaps;
};

#endif
