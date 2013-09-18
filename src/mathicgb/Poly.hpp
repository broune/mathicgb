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

  Poly(const PolyRing& ring): mRing(ring) {}

  void parse(std::istream &i); // reads into this, sorts terms
  void parseDoNotOrder(std::istream &i); // reads into this, does not sort terms
  void display(FILE* file, bool printComponent = true) const;
  void display(std::ostream& out, bool printComponent = true) const;
  void see(bool print_comp) const;

  class const_iterator {
  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef std::pair<const coefficient, const const_monomial> value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type* pointer; // todo: is this OK?
    typedef std::pair<const coefficient&, const const_monomial> reference;

    const_iterator() {}
    const_iterator& operator++() { ++ic; im += monsize; return *this; }

    coefficient coef() const {return *ic;}
    ConstMonoRef mono() const {return Monoid::toRef(&*im);}

    friend bool operator==(const const_iterator &a, const const_iterator &b);
    friend bool operator!=(const const_iterator &a, const const_iterator &b);

    NewConstTerm operator*() const {
      NewConstTerm t = {coef(), mono()};
      return t;
    }

  private:
    size_t monsize;
    std::vector<coefficient>::const_iterator ic;
    std::vector<exponent>::const_iterator im;
    friend class Poly;

    const_iterator(const Poly& f) : monsize(f.ring().maxMonomialSize()), ic(f.coeffs.begin()), im(f.monoms.begin()) {}
    const_iterator(const Poly& f,int) : ic(f.coeffs.end()), im() {}
  };

  const_iterator begin() const { return const_iterator(*this); }
  const_iterator end() const { return const_iterator(*this,1); }

  class MonoIterator {
    size_t monsize;
    std::vector<coefficient>::const_iterator ic;
    std::vector<exponent>::const_iterator im;
    friend class Poly;

    MonoIterator(const Poly& f) : monsize(f.ring().maxMonomialSize()), ic(f.coeffs.begin()), im(f.monoms.begin()) {}
    MonoIterator(const Poly& f,int) : ic(f.coeffs.end()), im(f.monoms.end()) {}

  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef ConstMonoRef value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef ConstMonoRef reference;

    MonoIterator() {}
    MonoIterator operator++() { ++ic; im += monsize; return *this; }
    bool operator==(const MonoIterator& it) const {return im == it.im;}
    bool operator!=(const MonoIterator& it) const {return im != it.im;}
    const value_type operator*() const {return Monoid::toRef(&*im);}
  };

  MonoIterator monoBegin() const { return MonoIterator(*this); }
  MonoIterator monoEnd() const { return MonoIterator(*this,1); }
  struct MonoRange {
    MonoIterator mBegin;
    const MonoIterator mEnd;
    MonoIterator begin() {return mBegin;}
    MonoIterator end() {return mEnd;}
  };
  MonoRange monoRange() const {
    MonoRange range = {monoBegin(), monoEnd()};
    return range;
  }

  /// Orders terms in descending order.
  void sortTermsDescending();

  /// Returns the coefficient of the given term.
  coefficient& coef(size_t index);

  /// Returns the coefficient of the given term.
  coefficient coef(size_t index) const;

  /// Returns the monomial of the given term.
  MonoRef mono(size_t index);

  /// Returns the monomial of the given term.
  ConstMonoRef mono(size_t index) const;

  /// Returns the coefficient of the leading term.
  coefficient leadCoef() const {return coef(0);}

  /// Returns the monomial of the leading term.
  ConstMonoRef leadMono() const {return mono(0);}

  /// Returns the monomial of the last term.
  const_monomial backMono() const;

  /// Appends the given term as the last term in the polynomial.
  void append(coefficient coef, ConstMonoRef mono);

  /// Hint that space for termCount terms is going to be needed. This serves
  /// the same purpose as std::vector<>::reserve.
  void reserve(size_t termCount);

  const coefficient* coefficientBegin() const {return coeffs.data();}

  void makeMonic();
  bool isMonic() const;

  bool isZero() const { return coeffs.empty(); }

  size_t termCount() const {return coeffs.size();}

  size_t getMemoryUse() const;

  void setToZero();

  Poly& operator=(const Poly& poly) {return *this = Poly(poly);}
  Poly& operator=(Poly&& poly);

  friend bool operator==(const Poly &a, const Poly &b);

  const PolyRing& ring() const {return mRing;}
  const Monoid& monoid() const {return ring().monoid();}

  bool termsAreInDescendingOrder() const;

private:
  const PolyRing& mRing;
  std::vector<coefficient> coeffs;
  std::vector<exponent> monoms;
};

std::ostream& operator<<(std::ostream& out, const Poly& p);

inline bool operator==(const Poly::const_iterator &a, const Poly::const_iterator &b)
{
  return a.ic == b.ic;
}
inline bool operator!=(const Poly::const_iterator &a, const Poly::const_iterator &b)
{
  return a.ic != b.ic;
}

inline void Poly::append(coefficient a, ConstMonoRef m) {
  coeffs.push_back(a);
  size_t len = ring().maxMonomialSize();
  auto& monoid = ring().monoid();
  const auto offset = monoms.size();
  monoms.resize(offset + monoid.entryCount());
  monoid.copy(m, *PolyRing::Monoid::MonoPtr(monoms.data() + offset));
}

MATHICGB_NAMESPACE_END
#endif
