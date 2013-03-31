#ifndef K_D_TREE_MODEL_GUARD
#define K_D_TREE_MODEL_GUARD

#include <string>
#include <vector>
#include <iostream>
#include <memtailor.h>
#include <mathic.h>

#include "PolyRing.hpp"

/** Helper class for KDTreeModel. */
template<bool UseDivMask, bool UseTreeDivMask, size_t LeafSize, bool AllowRemovals>
class KDTreeModelConfiguration;

template<bool UDM, bool UTDM, size_t LS, bool AR>
class MonTableKDTreeConfiguration {
 public:
  typedef int Exponent;
  typedef const_monomial Monomial;
  typedef Monomial Entry;

  MonTableKDTreeConfiguration
  (const PolyRing *R0,
     bool minimizeOnInsert,
     bool sortOnInsert,
     bool useDivisorCache,
     double rebuildRatio,
     size_t minRebuild):
        R(R0),
        _varCount(R->getNumVars()),
        _minimize_on_insert(minimizeOnInsert),
   _sortOnInsert(sortOnInsert),
   _useDivisorCache(useDivisorCache),
   _useAutomaticRebuild((rebuildRatio > 0.0 || minRebuild > 0) && UDM),
   _rebuildRatio(rebuildRatio),
   _minRebuild(minRebuild),
   _expQueryCount(0) {
     MATHICGB_ASSERT(rebuildRatio >= 0);
   }

  size_t getVarCount() const {return _varCount;}
  bool getSortOnInsert() const {return _sortOnInsert;}

  Exponent getExponent(const Monomial& m, size_t var) const {
    ++_expQueryCount;
    return R->monomialExponent(m, var);
  }

  bool divides(const Monomial& a, const Monomial& b) const {
    for (size_t var = 0; var < getVarCount(); ++var)
      if (getExponent(b, var) < getExponent(a, var))
        return false;
    return true;
  }

  bool isLessThan(const Monomial& a, const Monomial& b) const {
    for (size_t var = 0; var < getVarCount(); ++var) {
      if (getExponent(a, var) < getExponent(b, var))
        return true;
      if (getExponent(b, var) < getExponent(a, var))
        return false;
        }
    return false;
  }

  const PolyRing *getPolyRing() const {return R;}
  bool getUseDivisorCache() const {return _useDivisorCache;}
  bool getDoAutomaticRebuilds() const {return _useAutomaticRebuild;}
  double getRebuildRatio() const {return _rebuildRatio;}
  size_t getRebuildMin() const {return _minRebuild;}
  bool getMinimizeOnInsert() const {return _minimize_on_insert;}
  static const bool UseDivMask = UDM;
  static const bool UseTreeDivMask = UTDM;
  static const size_t LeafSize = LS;
  static const bool PackedTree = true;
  static const bool AllowRemovals = AR;

  unsigned long long getExpQueryCount() const {return _expQueryCount;}

 private:
  const PolyRing *R;
  const size_t _varCount;
  const bool _minimize_on_insert;
  const bool _sortOnInsert;
  const bool _useDivisorCache;
  const bool _useAutomaticRebuild;
  const double _rebuildRatio;
  const size_t _minRebuild;
  mutable unsigned long long _expQueryCount;
};

/** An instantiation of the capabilities of KDTree. */
template<bool UseDivMask, bool UseTreeDivMask, size_t LeafSize, bool AllowRemovals>
class MonTableKDTree {
 private:
  typedef MonTableKDTreeConfiguration<UseDivMask, UseTreeDivMask, LeafSize, AllowRemovals> C;
  typedef mic::KDTree<C> Finder;
 public:
  typedef typename Finder::Monomial Monomial;
  typedef typename Finder::Entry Entry;
  typedef C Configuration;

  MonTableKDTree(const PolyRing *R,
                 bool minimizeOnInsert,
                 bool sortOnInsert,
                 bool useDivisorCache,
                 double rebuildRatio,
                 size_t minRebuild):
        _finder(C(R, minimizeOnInsert, sortOnInsert, useDivisorCache, rebuildRatio, minRebuild))
  {
    MATHICGB_ASSERT(!UseTreeDivMask || UseDivMask);
  }

  MonTableKDTree(const Configuration &C) :
        _finder(C)
  {
    MATHICGB_ASSERT(!UseTreeDivMask || UseDivMask);
  }

  const C& getConfiguration() const {return _finder.getConfiguration();}
  C& getConfiguration() {return _finder.getConfiguration();}

  bool insert(const Entry& entry);
  template<class MultipleOutput>
  bool insert(const Entry& entry, MultipleOutput& removed);

  inline Entry* findDivisor(const Monomial& monomial) {
    return _finder.findDivisor(monomial);
  }
  inline const Entry* findDivisor(const Monomial& monomial) const {
    return _finder.findDivisor(monomial);
  }
  std::string getName() const;

  template<class DO>
  inline void findAllDivisors(const Monomial& monomial, DO& out) {
    _finder.findAllDivisors(monomial, out);
  }
  template<class DO>
  inline void findAllDivisors(const Monomial& monomial, DO& out) const {
    _finder.findAllDivisors(monomial, out);
  }

  unsigned long long getExpQueryCount() const {
    return _finder.getConfiguration().getExpQueryCount();
  }

  void getMonomials(std::vector<const_monomial>& monomials);

public:
  typedef size_t ValueType;

  const PolyRing * getPolyRing() const { return getConfiguration().getPolyRing(); }

  inline bool member(const_monomial t, ValueType & val) const {
    const Entry* entry = findDivisor(t);
    if (entry == 0)
      return false;
    val = 0;
    return true;
  }
  bool insert(const_monomial t, ValueType /* val */) { return insert(t); } // Only insert if not there

  size_t n_elems() const { return _finder.size(); }

  size_t getMemoryUse() const;

  void displayStats(std::ostream &o) const;

  struct Stats {
    size_t n_member;
    size_t n_inserts;  // includes koszuls
    size_t n_insert_already_there;
    size_t n_compares;
  };

  void getStats(Stats &stats) const { stats = stats_; }

 private:
  Finder _finder;
  Stats stats_;
};

class MonomialDeleter {
public:
  MonomialDeleter(const PolyRing& ring): mRing(ring) {}
  void push_back(const void* p) {
    mRing.freeMonomial(static_cast<exponent*>(const_cast<void*>(p)));
  }
  void push_back(ConstMonomial p) {
    mRing.freeMonomial(const_cast<exponent *>(p.unsafeGetRepresentation()));
  }
private:
  const PolyRing& mRing;
};

template<bool UDM, bool UTDM, size_t LS, bool AR>
inline bool MonTableKDTree<UDM, UTDM, LS, AR>::insert(const Entry& entry) {
  if (!getConfiguration().getMinimizeOnInsert()) {
    _finder.insert(entry);
    return true;
  }
  if (C::AllowRemovals) {
    if (findDivisor(entry) != 0)
      return false;
    MonomialDeleter mondelete(*getPolyRing());
    _finder.removeMultiples(entry, mondelete);
  } else {
    MATHICGB_ASSERT(findDivisor(entry) == 0);
    MonomialDeleter mondelete(*getPolyRing());
    // todo: can't do the below ASSERT as KD tree won't allow a
    // removal action when removals are not allowed. So there
    // should be a getMultiples() method that could be
    // used instead.
    //MATHICGB_ASSERT(!_finder.removeMultiples(entry, mondelete));
  }
  _finder.insert(entry);
  return true;
}

template<bool UDM, bool UTDM, size_t LS, bool AR>
template<class MultipleOutput>
inline bool MonTableKDTree<UDM, UTDM, LS, AR>::
insert(const Entry& entry, MultipleOutput& removed) {
  if (!getConfiguration().getMinimizeOnInsert()) {
    _finder.insert(entry);
    return true;
  }
  if (C::AllowRemovals) {
    if (findDivisor(entry) != _finder.end())
      return false;
    _finder.removeMultiples(entry, removed);
  } else {
    MATHICGB_ASSERT(findDivisor(entry) == _finder.end());
    // todo: can't do the below ASSERT as KD tree won't allow a
    // removal action when removals are not allowed. So there
    // should be a getMultiples() method that could be
    // used instead.
    //MATHICGB_ASSERT(!_finder.removeMultiples(entry, removed));
  }
  _finder.insert(entry);
  return true;
}

template<bool UDM, bool UTDM, size_t LS, bool AR>
inline std::string MonTableKDTree<UDM, UTDM, LS, AR>::getName() const {
  return _finder.getName() +
    (getConfiguration().getMinimizeOnInsert() ? " remin" : " nomin");
}

template<bool UDM, bool UTDM, size_t LS, bool AR>
size_t MonTableKDTree<UDM, UTDM, LS, AR>::getMemoryUse() const
{
  return _finder.getMemoryUse();
}

template <bool UDM, bool UTM, size_t LS, bool AR>
void MonTableKDTree<UDM,UTM,LS, AR>::displayStats(std::ostream &o) const
{
  o << "  #elements: " << n_elems() <<  std::endl;
  o << "  #calls member: " << stats_.n_member << std::endl;
  o << "  #calls insert, but already there: " << stats_.n_insert_already_there << std::endl;
  o << "  #compares: " << stats_.n_compares << std::endl;
}

namespace {
  class ForAllCopy {
  public:
    ForAllCopy(std::vector<const_monomial>& monomials):
      mMonomials(monomials) {}
    bool proceed(const_monomial m) {
      mMonomials.push_back(m);
      return true;
    }
  private:
    std::vector<const_monomial>& mMonomials;
  };
}

template<bool UDM, bool UTM, size_t LS, bool AR>
void MonTableKDTree<UDM, UTM, LS, AR>::getMonomials(std::vector<const_monomial>& monomials) {
  ForAllCopy copier(monomials);
  _finder.forAll(copier);
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
