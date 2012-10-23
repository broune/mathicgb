// Copyright 2011 Michael E. Stillman

#ifndef _poly_h_
#define _poly_h_

#include "PolyRing.hpp"
#include <vector>
#include <ostream>
#include <utility>
#include <cstdio>

class Poly {
public:
  Poly(const PolyRing& ring) : R(&ring) {MATHICGB_ASSERT(R != 0);}

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

    iterator(Poly& f) : monsize(f.getRing()->maxMonomialSize()), ic(f.coeffs.begin()), im(f.monoms.begin()) {}
    iterator(Poly& f,int) : ic(f.coeffs.end()), im() {}
  public:
    iterator() {}
    iterator operator++() { ++ic; im += monsize; return *this; }
    coefficient &getCoefficient() const { return *ic; }
    monomial getMonomial() const { return &*im; }
    size_t operator-(const iterator &b) const { return ic - b.ic; }
    friend bool operator==(const iterator &a, const iterator &b);
    friend bool operator!=(const iterator &a, const iterator &b);
    std::pair<coefficient&, monomial> operator*() const {
      return std::pair<coefficient&, monomial>(getCoefficient(), getMonomial());
    }
  };

  class const_iterator {
    // only for const objects...
    size_t monsize;
    std::vector<coefficient>::const_iterator ic;
    std::vector<exponent>::const_iterator im;
    friend class Poly;

    const_iterator(const Poly& f) : monsize(f.getRing()->maxMonomialSize()), ic(f.coeffs.begin()), im(f.monoms.begin()) {}
    const_iterator(const Poly& f,int) : ic(f.coeffs.end()), im() {}
  public:
    const_iterator() {}
    const_iterator operator++() { ++ic; im += monsize; return *this; }
    coefficient getCoefficient() const { return *ic; }
    const_monomial getMonomial() const { return &*im; }
    size_t operator-(const const_iterator &b) const { return ic - b.ic; }
    friend bool operator==(const const_iterator &a, const const_iterator &b);
    friend bool operator!=(const const_iterator &a, const const_iterator &b);
    std::pair<coefficient, const_monomial> operator*() const {
      return std::pair<coefficient, const_monomial>
        (getCoefficient(), getMonomial());
    }
  };

  void sortTermsDescending();

  void append(iterator &first, iterator &last);
  // Insert [first, last) onto the end of this.
  // This invalidates iterators on this.

  Poly *copy() const;

  // Cannot call this monomial() since that is already a type :-(
  monomial monomialAt(size_t index);
  const_monomial monomialAt(size_t index) const;
  coefficient& coefficientAt(size_t index);
  const coefficient coefficientAt(size_t index) const;

  /// all iterators are invalid after this
  void appendTerm(coefficient a, const_monomial m);

  /// Hint that space for termCount terms is going to be needed so the internal
  /// storage should be expanded to fit that many terms.
  void reserve(size_t spaceForThisManyTerms);

  const_iterator begin() const { return const_iterator(*this); }
  const_iterator end() const { return const_iterator(*this,1); }
  const_monomial backMonomial() const;

  iterator begin() { return iterator(*this); }
  iterator end() { return iterator(*this,1); }


  static Poly * add(const PolyRing *R,
                    iterator first1,
                    iterator last1,
                    iterator first2,
                    iterator last2,
                    size_t &n_compares);

  void multByTerm(coefficient a, const_monomial m);
  void multByMonomial(const_monomial m);
  void multByCoefficient(coefficient a);

  void makeMonic();
  bool isMonic() const;

  const_monomial getLeadMonomial() const { return &(monoms[0]); }
  const_coefficient getLeadCoefficient() const  { return coeffs[0]; }
  exponent getLeadComponent() const  { return R->monomialGetComponent(&(monoms[0])); }
  bool isZero() const { return coeffs.empty(); }
  size_t nTerms() const { return coeffs.size(); }

  size_t getMemoryUse() const;

  void setToZero();

  void copy(Poly &result) const;
  friend bool operator==(const Poly &a, const Poly &b);

  // deprecated
  const PolyRing *getRing() const { return R; }

  // use this instead of getRing()
  const PolyRing& ring() const {return *R;}

  void dump() const; // used for debugging

  bool termsAreInDescendingOrder() const;

private:
  const PolyRing *R;
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
  size_t len = R->monomialSize(m);
  exponent const * e = m.unsafeGetRepresentation();
  monoms.insert(monoms.end(), e, e + len);
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
