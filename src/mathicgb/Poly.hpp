// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_POLY_GUARD
#define MATHICGB_POLY_GUARD

#include "PolyRing.hpp"
#include "Range.hpp"
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

  Poly(const Poly& poly):
    mRing(poly.ring()), mCoefs(poly.mCoefs), mMonos(poly.mMonos)
  {}

  Poly(const Poly&& poly):
    mRing(poly.ring()),
    mCoefs(std::move(poly.mCoefs)),
    mMonos(std::move(poly.mMonos))
  {}

  const PolyRing& ring() const {return mRing;}
  const Monoid& monoid() const {return ring().monoid();}
  bool isZero() const {return mCoefs.empty();}
  size_t termCount() const {return mCoefs.size();}
  size_t getMemoryUse() const;

  /// Returns a polynomial whose terms have been permuted to be in
  /// descending order.
  ///
  /// Making the copy is not wasteful, because doing the permutation in-place
  /// would be require swapping monomials which is slow if they are large.
  /// The returned object is not copy (return value optimization) and using
  /// move assignment this code will only create the single copy of a Poly
  /// p that is necessary to avoid an in-place operation:
  ///
  ///   p = p.polyWithTermsDescending()
  Poly polyWithTermsDescending();

  /// Appends the given term as the last term in the polynomial.
  void append(const NewConstTerm& term) {
    MATHICGB_ASSERT(term.mono != nullptr);
    append(term.coef, *term.mono);
  }

  void append(coefficient coef, ConstMonoRef mono);

  /// Hint that space for the give number of terms is going to be needed.
  /// This serves the same purpose as std::vector<>::reserve.
  void reserve(size_t spaceForThisManyTerms) {
    mMonos.reserve(spaceForThisManyTerms * monoid().entryCount());
  }

  /// Makes the polynomial monic by multiplying by the multiplicative inverse
  /// of leadCoef(). Calling this method is an error if isZero().
  void makeMonic();

  void setToZero() {
    mCoefs.clear();
    mMonos.clear();
  }

  Poly& operator=(const Poly& poly) {return *this = Poly(poly);}
  Poly& operator=(Poly&& poly) {
    MATHICGB_ASSERT(&ring() == &poly.ring());
    mCoefs = std::move(poly.mCoefs);
    mMonos = std::move(poly.mMonos);
    return *this;
  }


  // *** Accessing the coefficients of the terms in the polynomial.

  /// Returns the coefficient of the given term.
  coefficient coef(size_t index) const {
    MATHICGB_ASSERT(index < termCount());
    return mCoefs[index];
  }

  /// Returns the coefficient of the leading term.
  coefficient leadCoef() const {
    MATHICGB_ASSERT(!isZero());
    return mCoefs.front();
  }

  /// Returns true if the polynomial is monic. A polynomial is monic if
  /// the coefficient of the leading monomial is 1. If you are asking this
  /// question about a polynomial, that likely means that you are expecting
  /// the polynomial not to be zero. So it is an error to ask if the zero
  /// polynomial is monic - you'll get an assert to help pinpoint the error.
  bool isMonic() const {
    MATHICGB_ASSERT(!isZero());
    return ring().coefficientIsOne(leadCoef());
  }

  typedef std::vector<coefficient>::const_iterator ConstCoefIterator;
  typedef Range<ConstCoefIterator> ConstCoefIteratorRange;

  ConstCoefIterator coefBegin() const {return mCoefs.begin();}
  ConstCoefIterator coefEnd() const {return mCoefs.end();}
  ConstCoefIteratorRange coefRange() const {
    return range(coefBegin(), coefEnd());
  }


  // *** Accessing the monomials of the terms in the polynomial

  /// Returns the monomial of the given term.
  ConstMonoRef mono(size_t index) const {
    MATHICGB_ASSERT(index < termCount());
    return Monoid::toRef(&mMonos[index * monoid().entryCount()]);
  }

  /// Returns the monomial of the leading term.
  ConstMonoRef leadMono() const {
    MATHICGB_ASSERT(!isZero());
    return Monoid::toRef(&mMonos.front());
  }

  /// Returns the monomial of the last term.
  ConstMonoRef backMono() const {
    MATHICGB_ASSERT(!isZero());
    return Monoid::toRef(&mMonos[(termCount() - 1) * monoid().entryCount()]);
  }

  /// Returns true if the terms are in descending order. The terms are in
  /// descending order when mono(0) >= mono(1) >= ... >= backMono.
  bool termsAreInDescendingOrder() const;

  class ConstMonoIterator {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef ConstMonoRef value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef ConstMonoRef reference;

    ConstMonoIterator() {}

    ConstMonoIterator& operator++() {
      mIt += mEntryCount;
      return *this;
    }

    ConstMonoRef operator*() const {return Monoid::toRef(&*mIt);}

    bool operator==(const ConstMonoIterator& it) const {return mIt == it.mIt;}
    bool operator!=(const ConstMonoIterator& it) const {return mIt != it.mIt;}

  private:
    friend class Poly;
    typedef std::vector<exponent>::const_iterator Iterator;

    ConstMonoIterator(const Monoid& monoid, Iterator it):
      mEntryCount(monoid.entryCount()),
      mIt(it)
    {}

    size_t mEntryCount;
    Iterator mIt;
  };

  typedef Range<ConstMonoIterator> ConstMonoIteratorRange;

  ConstMonoIterator monoBegin() const {
    return ConstMonoIterator(monoid(), mMonos.begin());
  }

  ConstMonoIterator monoEnd() const {
    return ConstMonoIterator(monoid(), mMonos.end());
  }

  ConstMonoIteratorRange monoRange() const {
    return range(monoBegin(), monoEnd());
  }


  // *** Iteration through terms

  class ConstTermIterator {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef NewConstTerm value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    ConstTermIterator() {}

    ConstTermIterator& operator++() {
      ++mIt;
      return *this;
    }

    value_type operator*() const {
      auto pair = *mIt;
      NewConstTerm term = {pair.first, pair.second};
      return term;
    }

    bool operator==(const ConstTermIterator& it) const {return mIt == it.mIt;}
    bool operator!=(const ConstTermIterator& it) const {return mIt != it.mIt;}

    coefficient coef() const {return (*mIt).first;}
    ConstMonoRef mono() const {return (*mIt).second;}

  private:
    friend class Poly;
    typedef Zip<ConstCoefIterator, ConstMonoIterator> Iterator;
    ConstTermIterator(const Iterator& it): mIt(it) {}

    Iterator mIt;
  };

  NewConstTerm term(size_t index) const {
    NewConstTerm t = {coef(index), mono(index).ptr()};
    return t;
  }

  typedef Range<ConstTermIterator> ConstTermIteratorRange;

  ConstTermIterator begin() const {return makeZip(coefBegin(), monoBegin());}
  ConstTermIterator end() const {return makeZip(coefEnd(), monoEnd());}
  ConstTermIteratorRange termRange() const {return range(begin(), end());}

private:
  friend bool operator==(const Poly &a, const Poly &b);

  const PolyRing& mRing;
  std::vector<coefficient> mCoefs;
  std::vector<exponent> mMonos;
};

// This is inline since it is performance-critical.
inline void Poly::append(coefficient a, ConstMonoRef m) {
  mCoefs.push_back(a);

  const auto offset = mMonos.size();
  mMonos.resize(offset + monoid().entryCount());
  monoid().copy(m, *PolyRing::Monoid::MonoPtr(mMonos.data() + offset));
}

MATHICGB_NAMESPACE_END
#endif
