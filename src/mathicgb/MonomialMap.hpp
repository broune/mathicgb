#ifndef MATHICGB_MONOMIAL_MAP_GUARD
#define MATHICGB_MONOMIAL_MAP_GUARD

#include "PolyRing.hpp"
#include <limits>

#if defined(MATHICGB_USE_STD_HASH) && defined(MATHICGB_USE_CUSTOM_HASH)
#error Only select one kind of hash table
#endif

// set default hash table type if nothing has been specified
#if !defined(MATHICGB_USE_STD_HASH) && !defined(MATHICGB_USE_CUSTOM_HASH)
#define MATHICGB_USE_CUSTOM_HASH
#endif

#ifdef MATHICGB_USE_CUSTOM_HASH
#include <memtailor.h>
#include <vector>
#include <algorithm>
#elif defined(MATHICGB_USE_STD_HASH)
#include <unordered_map>
#else
#include <map>
#endif

/** A mapping from monomials to MapToType. The interface is the same as
  for STL maps. The interface is not complete. Add methods you need when
  you need it.
*/
template<class MapToType>
class MonomialMap;

namespace MonomialMapInternal {
  // The map type MapClass is defined here. This is here in
  // this namespace to avoid cluttering the class definition with
  // a lot of ifdef's.

#ifdef MATHICGB_USE_CUSTOM_HASH
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
    MapClass(const PolyRing& ring):
      mRing(ring),
      mTable(),
      mNodeAlloc(sizeof(Node) - sizeof(exponent) + ring.maxMonomialByteSize())
    {
      mGrowWhenThisManyEntries = 0;
      mTableSize = mTable.size();
      mEntryCount = 0;
      growTable();
    }

    Map& map() {return *this;}
    const Map& map() const {return *this;}
    const PolyRing& ring() const {return mRing;}

    mapped_type* find(const_monomial mono) {
      const HashValue monoHash = mRing.monomialHashValue(mono);
      Node* node = entry(hashToIndex(monoHash));
      for (; node != 0; node = node->next) {
        if (monoHash == mRing.monomialHashValue(node->mono) &&
          mRing.monomialEqualHintTrue(mono, node->mono)
        ) {
          return &node->value;
        }
      }
      return 0;
    }

    mapped_type* findProduct(const_monomial a, const_monomial b) {
      const HashValue abHash = mRing.monomialHashOfProduct(a, b);
      Node* node = entry(hashToIndex(abHash));
      for (; node != 0; node = node->next) {
        if (abHash == mRing.monomialHashValue(node->mono) &&
          mRing.monomialIsProductOfHintTrue(a, b, node->mono)
        ) {
          return &node->value;
        }
      }
      return 0;
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

#elif defined(MATHICGB_USE_STD_HASH)
  class Hash {
  public:
    Hash(const PolyRing& ring): mRing(ring) {}
    size_t operator()(const_monomial m) const {
      return mRing.monomialHashValue(m);
    }
    const PolyRing& ring() const {return mRing;}

  private:
    const PolyRing& mRing;
  };

  class Equal {
  public:
    Equal(const PolyRing& ring): mRing(ring) {}
    size_t operator()(const_monomial a, const_monomial b) const {
      return mRing.monomialEQ(a, b);
    }
    const PolyRing& ring() const {return mRing;}

  private:
    const PolyRing& mRing;
  };

  // This has to be a template to support rebind().
  template<class T>
  class TemplateHashAllocator {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
    typedef const T* const_pointer;
    typedef const T& const_reference;
    typedef ::size_t size_type;
    typedef ::ptrdiff_t difference_type;

    template<class T2>
    struct rebind {
      typedef TemplateHashAllocator<T2> other;
    };

    TemplateHashAllocator() {}
    TemplateHashAllocator(memt::Arena& arena): mArena(arena) {}
    template<class X>
    TemplateHashAllocator(const TemplateHashAllocator<X>& a):
      mArena(a.arena()) {}
    TemplateHashAllocator(const TemplateHashAllocator<T>& a):
      mArena(a.arena()) {}

    pointer address(reference x) {return &x;}
    const_pointer address(const_reference x) const {return &x;}
    pointer allocate(size_type n, void* hint = 0) {
      return static_cast<pointer>(mArena.alloc(sizeof(T) * n));
    }
    void deallocate(pointer p, size_t n) {}
    size_type max_size() const {return std::numeric_limits<size_type>::max();}
    void construct(pointer p, const_reference val) {new (p) T(val);}
    void destroy(pointer p) {p->~T();}
    memt::Arena& arena() const {return mArena;}

  private:
    mutable memt::Arena& mArena;
  };

  template<class MTT>
  class MapClass {
  public:
    typedef TemplateHashAllocator<std::pair<const const_monomial, MTT>>
      HashAllocator;
    typedef std::unordered_map<const_monomial, MTT, Hash, Equal, HashAllocator>
      Map;
    typedef typename Map::iterator iterator;
    typedef typename Map::const_iterator const_iterator;
    typedef typename Map::value_type value_type;

    MapClass(const PolyRing& ring):
      mArena(),
      mMap(100, Hash(ring), Equal(ring), HashAllocator(mArena)),
      mTmp(ring.allocMonomial())
    {
      mMap.max_load_factor(0.3f);
    }

    ~MapClass() {
      ring().freeMonomial(mTmp);
    }

    Map& map() {return mMap;}
    const Map& map() const {return mMap;}
    const PolyRing& ring() const {return mMap.key_eq().ring();}

    value_type* findProduct(const_monomial a, const_monomial b) {
      ring().setGivenHash(mTmp, ring().monomialHashOfProduct(a, b));
      size_t bucket = mMap.bucket(mTmp);
      auto stop = mMap.end(bucket);
      for (auto it = mMap.begin(bucket); it != stop; ++it)
        if (ring().monomialIsProductOf(a, b, it->first))
          return &*it;
      return 0;
    }

  private:
    memt::Arena mArena;
    Map mMap;
    monomial mTmp;
  };
#elif defined(MATHICGB_USE_CUSTOM_HASH)
#else
  /// We need SOME ordering to make std::map work.
  class ArbitraryOrdering {
  public:
    ArbitraryOrdering(const PolyRing& ring): mRing(ring) {}
    bool operator()(const_monomial a, const_monomial b) const {
      return mRing.monomialLT(a, b);
    }
    const PolyRing& ring() const {return mRing;}

  private:
    const PolyRing& mRing;
  };

  template<class MTT>
  class MapClass {
  public:
    typedef std::map<const_monomial, MTT, ArbitraryOrdering> Map;

    MapClass(const PolyRing& ring): mMap(ArbitraryOrdering(ring)) {}
    Map& map() {return mMap;}
    const Map& map() const {return mMap;}
    const PolyRing& ring() const {return mMap.key_cmp().ring();}

  private:
    Map mMap;
  };
#endif
}

template<class MTT>
class MonomialMap {
public:
  typedef MTT mapped_type;
  typedef MonomialMapInternal::MapClass<MTT> MapType;

  //typedef typename MapType::Map::iterator iterator;
  //typedef typename MapType::Map::const_iterator const_iterator;
  typedef typename MapType::Map::value_type value_type;

  MonomialMap(const PolyRing& ring): mMap(ring) {}

/*  iterator begin() {return map().begin();}
  const_iterator begin() const {return map().begin();}
  iterator end() {return map().end();}
  const_iterator end() const {return map().end();}*/

  mapped_type* find(const_monomial m) {return mMap.find(m);}

  mapped_type* findProduct(const_monomial a, const_monomial b) {
    return mMap.findProduct(a, b);
  }

  const mapped_type* find(const_monomial m) const {
    return const_cast<MonomialMap&>(*this).find(m);
  }

  const mapped_type* findProduct(const_monomial a, const_monomial b) const {
    return const_cast<MonomialMap&>(*this).findProduct(a, b);
  }

  size_t size() const {return map().size();}
  void clear() {map().clear();}
  void insert(const value_type& val) {
    map().insert(val);
  }

  inline const PolyRing& ring() const {return mMap.ring();}

private:
  typename MapType::Map& map() {return mMap.map();}
  const typename MapType::Map& map() const {return mMap.map();}

  MapType mMap;
};

#endif
