#ifndef DIV_ARRAY_MODEL_GUARD
#define DIV_ARRAY_MODEL_GUARD

#include <string>
#include <vector>
#include <iostream>

#include <mathic.h>
#include "PolyRing.hpp"

/** Helper class for MonTableDivList. */
template<bool UseLinkedList, bool UseDivMask>
class MonTableDivListConfiguration;

template<bool ULL, bool UDM>
class MonTableDivListConfiguration {
public:
  typedef int Exponent;
  typedef const_monomial Monomial;
  typedef Monomial Entry;

  MonTableDivListConfiguration
  (const PolyRing *R0,
   bool minimizeOnInsert,
   bool moveDivisorToFront,
   bool sortOnInsert,
   double rebuildRatio,
   size_t minRebuild):
    R(R0),
    _varCount(R->getNumVars()),
        _minimizeOnInsert(minimizeOnInsert),
        _moveDivisorToFront(moveDivisorToFront),
  _sortOnInsert(sortOnInsert),
  _useAutomaticRebuild((rebuildRatio > 0.0 || minRebuild > 0) && UDM),
  _rebuildRatio(rebuildRatio),
  _minRebuild(minRebuild),
  _expQueryCount(0) {}

  static const bool UseLinkedList = ULL;
  static const bool UseDivMask = UDM;

  const PolyRing *getPolyRing() const {return R;}
  bool getDoAutomaticRebuilds() const {return _useAutomaticRebuild;}
  double getRebuildRatio() const {return _rebuildRatio;}
  size_t getRebuildMin() const {return _minRebuild;}
  bool getSortOnInsert() const {return _sortOnInsert;}
  bool getMinimizeOnInsert() const {return _minimizeOnInsert;}
  bool getMoveDivisorToFront() const {return _moveDivisorToFront;}

  size_t getVarCount() const {return _varCount;}

  NO_PINLINE Exponent getExponent(const Monomial& m, size_t var) const {
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

  unsigned long long getExpQueryCount() const {return _expQueryCount;}

 private:
  const PolyRing *R;
  const size_t _varCount;
  const bool _minimizeOnInsert;
  const bool _moveDivisorToFront;
  const bool _sortOnInsert;
  const bool _useAutomaticRebuild;
  const double _rebuildRatio;
  const size_t _minRebuild;
  mutable unsigned long long _expQueryCount;
};

template<bool UseLinkedList, bool UseDivMask>
class MonTableDivList;

/** An instantiation of the capabilities of DivList. */
template<bool ULL, bool UDM>
class MonTableDivList {
 private:
  typedef MonTableDivListConfiguration<ULL, UDM> C;
  typedef mic::DivList<C> Finder;
 public:
  typedef typename Finder::iterator iterator;
  typedef typename Finder::const_iterator const_iterator;
  typedef typename Finder::Monomial Monomial;
  typedef typename Finder::Entry Entry;
  typedef C Configuration;

  MonTableDivList(const PolyRing *R,
                          bool minimizeOnInsert,
                          bool moveDivisorToFront,
                          bool sortOnInsert,
              double rebuildRatio,
              size_t minRebuild):
        _finder(C(R, minimizeOnInsert, moveDivisorToFront, sortOnInsert, rebuildRatio, minRebuild))
  {
    ASSERT(!sortOnInsert || !moveDivisorToFront);
  }

  MonTableDivList(const Configuration &C) :
        _finder(C)
  {
    ASSERT(!C.getSortOnInsert() || !C.getMoveDivisorToFront());
  }

  const C& getConfiguration() const {return _finder.getConfiguration();}
  C& getConfiguration() {return _finder.getConfiguration();}

  bool getMinimizeOnInsert() const {return getConfiguration().getMinimizeOnInsert();}
  bool getMoveDivisorToFront() const {return getConfiguration().getMoveDivisorToFront();}

  bool insert(const Entry& entry);
  template<class MultipleOutput>
  bool insert(const Entry& entry, MultipleOutput& removed);

  Entry* findDivisor(const Monomial& monomial) {
    return _finder.findDivisor(monomial);
  }

  const Entry* findDivisor(const Monomial& monomial) const {
    return const_cast<MonTableDivList<ULL, UDM>&>(*this).findDivisor(monomial);
  }

  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out) {
    _finder.findAllDivisors(monomial, out);
  }
  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out) const {
    _finder.findAllDivisors(monomial, out);
  }

  std::string getName() const;


  iterator begin() {return _finder.begin();}
  const_iterator begin() const {return _finder.begin();}
  iterator end() {return _finder.end();}
  const_iterator end() const {return _finder.end();}
  size_t size() const {return _finder.size();}

  unsigned long long getExpQueryCount() const {
    return _finder.getConfiguration().getExpQueryCount();
  }

  void getMonomials(std::vector<const_monomial>& monomials);

public:
  typedef size_t ValueType;

  const PolyRing * getPolyRing() const { return getConfiguration().getPolyRing(); }

  bool member(const_monomial t, ValueType & val) const {
      const Entry* i = findDivisor(t);
      if (i == 0) return false;
      val = 0;
      return true;
  }
  bool insert(const_monomial t, ValueType /* val */) { return insert(t); } // Only insert if not there

  void display(std::ostream &o, int level) const;

  void dump(int level) const;

  size_t n_elems() const { return size(); }

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

template<bool ULL, bool UDM>
inline bool MonTableDivList<ULL, UDM>::insert(const Entry& entry) {
  if (!getMinimizeOnInsert()) {
    _finder.insert(entry);
    return true;
  }
  if (findDivisor(entry) != 0)
    return false;
  bool hasMultiples = _finder.removeMultiples(entry);
  _finder.insert(entry);
  if (getMoveDivisorToFront() && hasMultiples) {
    iterator it = _finder.end();
    _finder.moveToFront(--it);
  }
  return true;
}

template<bool ULL, bool UDM>
template<class MO>
inline bool MonTableDivList<ULL, UDM>::insert(const Entry& entry, MO& out) {
  if (!getMinimizeOnInsert()) {
    _finder.insert(entry);
    return true;
  }
  if (findDivisor(entry) != _finder.end())
    return false;
  bool hasMultiples = _finder.removeMultiples(entry, out);
  _finder.insert(entry);
  if (getMoveDivisorToFront() && hasMultiples) {
    iterator it = _finder.end();
    _finder.moveToFront(--it);
  }
  return true;
}

template<bool ULL, bool UDM>
inline std::string MonTableDivList<ULL, UDM>::getName() const {
  return _finder.getName() +
    (getMinimizeOnInsert() ? " remin" : " nomin") +
    (getMoveDivisorToFront() ? " toFront" : "");
}

template <bool UDL, bool UDM>
void MonTableDivList<UDL,UDM>::display(std::ostream &o, int /* level */) const
{
  const_iterator e = _finder.end();

  for (const_iterator i=_finder.begin(); i != e; ++i)
    {
      getPolyRing()->monomialDisplay(o, *i, false);
      o << "  ";
    }
  o << std::endl;
}

template <bool UDL, bool UDM>
void MonTableDivList<UDL,UDM>::dump(int level) const
{
  displayStats(std::cerr);
  if (level > 0) display(std::cerr, level-1);
}

template <bool UDL, bool UDM>
size_t MonTableDivList<UDL,UDM>::getMemoryUse() const
{
  return _finder.getMemoryUse();
}

template <bool UDL, bool UDM>
void MonTableDivList<UDL,UDM>::displayStats(std::ostream &o) const
{
  o << "  #elements: " << n_elems() <<  std::endl;
  o << "  #calls member: " << stats_.n_member << std::endl;
  o << "  #calls insert, but already there: " << stats_.n_insert_already_there << std::endl;
  o << "  #compares: " << stats_.n_compares << std::endl;
}

template<bool UDL, bool UDM>
void MonTableDivList<UDL, UDM>::getMonomials(std::vector<const_monomial>& monomials) {
  monomials.insert(monomials.begin(), begin(), end());
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
