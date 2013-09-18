// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_POLY_GUARD
#define MATHICGB_POLY_GUARD

#include "PolyRing.hpp"
#include <vector>
#include <ostream>
#include <utility>
#include <cstdio>
#include <iterator>

MATHICGB_NAMESPACE_BEGIN

class Poly {
public:
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::Mono Mono;
  typedef Monoid::MonoRef MonoRef;
  typedef Monoid::ConstMonoRef ConstMonoRef;
  typedef Monoid::MonoPtr MonoPtr;
  typedef Monoid::ConstMonoPtr ConstMonoPtr;

  Poly(const PolyRing& ring) : mRing(ring) {}

  void parse(std::istream &i); // reads into this, sorts terms
  void parseDoNotOrder(std::istream &i); // reads into this, does not sort terms
  void display(FILE* file, bool printComponent = true) const;
  void display(std::ostream& out, bool printComponent = true) const;
  void see(bool print_comp) const;

  class iterator {
    // only for const objects...
    size_t monsize;
    std::vector<coefficient>::iterator ic;
    std::vector<exponent>::iterator im;
    friend class Poly;

    iterator(Poly& f) : monsize(f.ring().maxMonomialSize()), ic(f.coeffs.begin()), im(f.monoms.begin()) {}
    iterator(Poly& f,int) : ic(f.coeffs.end()), im() {}
  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef std::pair<coefficient, const const_monomial> value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type* pointer; // todo: is this OK?
    typedef std::pair<coefficient&, const const_monomial> reference;

    iterator() {}
    iterator operator++() { ++ic; im += monsize; return *this; }
    coefficient &getCoefficient() const { return *ic; }
    monomial getMonomial() const { return &*im; }
    MonoRef mono() const {return getMonomial();}
    size_t operator-(const iterator &b) const { return ic - b.ic; }
    friend bool operator==(const iterator &a, const iterator &b);
    friend bool operator!=(const iterator &a, const iterator &b);
    reference operator*() const {
      return std::pair<coefficient&, monomial>(getCoefficient(), getMonomial());
    }
    const_term term() const {
      const_term t;
      t.monom = getMonomial();
      t.coeff = getCoefficient();
      return t;
    }
  };

  class const_iterator {
    // only for const objects...
    size_t monsize;
    std::vector<coefficient>::const_iterator ic;
    std::vector<exponent>::const_iterator im;
    friend class Poly;

    const_iterator(const Poly& f) : monsize(f.ring().maxMonomialSize()), ic(f.coeffs.begin()), im(f.monoms.begin()) {}
    const_iterator(const Poly& f,int) : ic(f.coeffs.end()), im() {}
  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef std::pair<const coefficient, const const_monomial> value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type* pointer; // todo: is this OK?
    typedef std::pair<const coefficient&, const const_monomial> reference;

    const_iterator() {}
    const_iterator operator++() { ++ic; im += monsize; return *this; }
    coefficient getCoefficient() const { return *ic; }
    const_monomial getMonomial() const { return &*im; }
    ConstMonoRef mono() const {return getMonomial();}
    size_t operator-(const const_iterator &b) const { return ic - b.ic; }
    friend bool operator==(const const_iterator &a, const const_iterator &b);
    friend bool operator!=(const const_iterator &a, const const_iterator &b);
    const value_type operator*() const {
      return std::pair<coefficient, const_monomial>
        (getCoefficient(), getMonomial());
    }
    NewConstTerm term() const {
      NewConstTerm t;
      t.mono = getMonomial();
      t.coef = getCoefficient();
      return t;
    }
  };

  void sortTermsDescending();

  void append(iterator &first, iterator &last);
  // Insert [first, last) onto the end of this.
  // This invalidates iterators on this.

  //Poly *copy() const;

  // Cannot call this monomial() since that is already a type :-(
  monomial monomialAt(size_t index);
  const_monomial monomialAt(size_t index) const;
  coefficient& coefficientAt(size_t index);
  const coefficient coefficientAt(size_t index) const;

  /// all iterators are invalid after this
  void appendTerm(coefficient a, const_monomial m);
  void appendTerm(coefficient a, PolyRing::Monoid::ConstMonoRef m);

  /// Hint that space for termCount terms is going to be needed so the internal
  /// storage should be expanded to fit that many terms.
  void reserve(size_t spaceForThisManyTerms);

  const_iterator begin() const { return const_iterator(*this); }
  const_iterator end() const { return const_iterator(*this,1); }
  const_monomial backMonomial() const;

  iterator begin() { return iterator(*this); }
  iterator end() { return iterator(*this,1); }

  const coefficient* coefficientBegin() const {return coeffs.data();}

  void multByTerm(coefficient a, const_monomial m);
  void multByMonomial(const_monomial m);
  void multByCoefficient(coefficient a);

  void makeMonic();
  bool isMonic() const;

  const_monomial getLeadMonomial() const { return &(monoms[0]); }
  const_coefficient getLeadCoefficient() const  { return coeffs[0]; }
  exponent getLeadComponent() const  { return ring().monomialGetComponent(&(monoms[0])); }
  bool isZero() const { return coeffs.empty(); }

  //size_t nTerms() const { return coeffs.size(); } /// @todo: deprecated
  size_t termCount() const {return coeffs.size();}

  size_t getMemoryUse() const;

  void setToZero();

  //void copy(Poly &result) const;
  Poly& operator=(const Poly& poly) {return *this = Poly(poly);}
  Poly& operator=(Poly&& poly);

  friend bool operator==(const Poly &a, const Poly &b);

  const PolyRing& ring() const {return mRing;}
  const Monoid& monoid() const {return ring().monoid();}

  //void dump() const; // used for debugging

  bool termsAreInDescendingOrder() const;

private:
  const PolyRing& mRing;
  std::vector<coefficient> coeffs;
  std::vector<exponent> monoms;
};

std::ostream& operator<<(std::ostream& out, const Poly& p);

inline bool operator==(const Poly::iterator &a, const Poly::iterator &b)
{
  return a.ic == b.ic;
}
inline bool operator!=(const Poly::iterator &a, const Poly::iterator &b)
{
  return a.ic != b.ic;
}

inline bool operator==(const Poly::const_iterator &a, const Poly::const_iterator &b)
{
  return a.ic == b.ic;
}
inline bool operator!=(const Poly::const_iterator &a, const Poly::const_iterator &b)
{
  return a.ic != b.ic;
}

inline void Poly::appendTerm(coefficient a, const_monomial m)
{
  // the monomial will be copied on.
  coeffs.push_back(a);
  size_t len = ring().maxMonomialSize();
  exponent const * e = m.unsafeGetRepresentation();
  monoms.insert(monoms.end(), e, e + len);
}

inline void Poly::appendTerm(coefficient a, PolyRing::Monoid::ConstMonoRef m) {
  coeffs.push_back(a);
  size_t len = ring().maxMonomialSize();
  auto& monoid = ring().monoid();
  const auto offset = monoms.size();
  monoms.resize(offset + monoid.entryCount());
  monoid.copy(m, *PolyRing::Monoid::MonoPtr(monoms.data() + offset));
}

MATHICGB_NAMESPACE_END
#endif
