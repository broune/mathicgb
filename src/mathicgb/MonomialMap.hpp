#ifndef MATHICGB_MONOMIAL_MAP_GUARD
#define MATHICGB_MONOMIAL_MAP_GUARD

#include "PolyRing.hpp"
#include <limits>

#define MATHICGB_USE_STD_HASH

#ifdef MATHICGB_USE_STD_HASH
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

#ifndef MATHICGB_USE_STD_HASH
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

#else
  class Hash {
  public:
    Hash(const PolyRing& ring): mRing(ring) {}
    size_t operator()(const_monomial m) const {
      MATHICGB_ASSERT(mRing.hashValid(m));
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

    MapClass(const PolyRing& ring):
      mArena(),
      mMap(100, Hash(ring), Equal(ring), HashAllocator(mArena))
    {
      mMap.max_load_factor(0.3f);
    }

    Map& map() {return mMap;}
    const Map& map() const {return mMap;}
    const PolyRing& ring() const {return mMap.key_eq().ring();}

  private:
    memt::Arena mArena;
    Map mMap;
  };
#endif
}

template<class MTT>
class MonomialMap {
public:
  typedef MTT mapped_type;
  typedef MonomialMapInternal::MapClass<MTT> MapType;

  typedef typename MapType::Map::iterator iterator;
  typedef typename MapType::Map::const_iterator const_iterator;
  typedef typename MapType::Map::value_type value_type;

  MonomialMap(const PolyRing& ring): mMap(ring) {}

  iterator begin() {return map().begin();}
  const_iterator begin() const {return map().begin();}
  iterator end() {return map().end();}
  const_iterator end() const {return map().end();}
  iterator find(const_monomial m) {return map().find(m);}
  const_iterator find(const_monomial m) const {return map().find(m);}

  size_t size() const {return map().size();}
  void clear() {map().clear();}
  std::pair<iterator, bool> insert(const value_type& val) {
    return map().insert(val);
  }

  inline const PolyRing& ring() const {return mMap.ring();}

private:
  typename MapType::Map& map() {return mMap.map();}
  const typename MapType::Map& map() const {return mMap.map();}

  MapType mMap;
};



#endif
