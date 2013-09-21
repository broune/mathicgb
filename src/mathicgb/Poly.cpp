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
  auto greaterOrEqual = [&](const NewConstTerm& a, const NewConstTerm& b) {
    return monoid().lessThan(*b.mono, *a.mono);
  };
  auto terms = rangeToVector(*this);
  std::sort(std::begin(terms), std::end(terms), greaterOrEqual);

  // *** Make a new polynomial with terms in that order
  Poly poly(ring());
  poly.append(terms, termCount());

  MATHICGB_ASSERT(poly.termsAreInDescendingOrder());

  // This return statements causes no copy. The return value optimization
  // will be used at the option of the compiler. If a crappy compiler gets
  // that wrong, poly will be treated as an r-value, which is to say that
  // this code becomes equivalent to return std::move(poly). That happens
  // because poly is a local variable being returned, so the standard
  // allows movement out of poly in this particular situation - that is
  // safe/reasonable because the very next thing that will happen to poly is
  // that it will get destructed, so anyone in a position to know that the
  // contents of poly had been moved out would then also be using
  // a pointer to the now-invalid poly object, which invokes undefined
  // behavior anyway.
  //
  // Capturing the returned Poly into another poly also will not cause a copy.
  // Consider the code
  //
  //   Poly p1(q.polyWithTermsDescending());
  //   Poly p2(ring);
  //   p2 = q.polyWithTermsDescending();
  //   
  // The return value is an unnamed temporary that is constructed/assigned
  // into p1/p2. Unnamed temporaries are r-values, so the guts of those
  // temporary objects will be moved into p1/p2. There will be no copy.
  //
  // (Of course there is the one copy that we are doing further up into poly
  // to avoid doing an in-place sort. The point is that there are no further
  // copies than that one.)
  return poly;
}

void Poly::makeMonic() {
  MATHICGB_ASSERT(!isZero());
  if (isMonic())
    return;
  auto multiply = leadCoef();
  ring().coefficientReciprocalTo(multiply);
  for (auto& coef : mCoefs)
    ring().coefficientMultTo(coef, multiply);
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
