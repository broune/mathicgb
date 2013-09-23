// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "Poly.hpp"

#include <algorithm>

MATHICGB_NAMESPACE_BEGIN

Poly Poly::polyWithTermsDescending() {
  // *** Sort terms in descending order of monomial.
  // It would be possible but cumbersome to implement a sort directly
  // on mMonos. That way no allocation would need to happen, however
  // it is not clear that that would be any faster, since swapping around
  // monomials in-place is slow. Swapping terms is faster, since terms
  // just refer to the monomials. This way is also easier to implement.
  //
  /// @todo: make a separate TermSorter object that allows the terms vector
  /// to be reused between sorts. This should matter for sorting input
  /// ideals where there might be a lot of polynomials to go through.
  /// That way terms would not need to be allocated anew for each polynomial.
  auto greaterOrEqual = [&](const NewConstTerm& a, const NewConstTerm& b) {
    return monoid().lessThan(*b.mono, *a.mono);
  };
  auto terms = rangeToVector(*this);
  std::sort(std::begin(terms), std::end(terms), greaterOrEqual);

  // *** Make a new polynomial with terms in that order
  Poly poly(ring());
  poly.reserve(termCount());
  poly.append(terms);

  MATHICGB_ASSERT(poly.termsAreInDescendingOrder());
  MATHICGB_ASSERT(poly.termCount() == termCount());
  return poly;
}

void Poly::makeMonic() {
  MATHICGB_ASSERT(!isZero());
  if (isMonic())
    return;
  auto multiplier = field().inverse(leadCoef());
  for (auto& coef : mCoefs)
    coef = field().product(coef, multiplier);
  MATHICGB_ASSERT(isMonic());
}

size_t Poly::getMemoryUse() const {
  return 
    sizeof(mCoefs.front()) * mCoefs.capacity() +
    sizeof(mMonos.front()) * mMonos.capacity();
}

bool Poly::termsAreInDescendingOrder() const {
  auto greaterThanOrEqual = [&](ConstMonoRef a, ConstMonoRef b) {
    return !monoid().lessThan(a, b);
  };
  return std::is_sorted(monoBegin(), monoEnd(), greaterThanOrEqual);
}

bool operator==(const Poly& a, const Poly& b) {
  MATHICGB_ASSERT(&a.ring() == &b.ring());
  if (a.termCount() != b.termCount())
    return false;

  for (const auto& coef : zip(a.coefRange(), b.coefRange()))
    if (coef.first != coef.second)
      return false;

  const auto& monoid = a.ring().monoid();
  for (const auto& mono : zip(a.monoRange(), b.monoRange()))
    if (!monoid.equal(mono.first, mono.second))
      return false;

  return true;
}

MATHICGB_NAMESPACE_END
