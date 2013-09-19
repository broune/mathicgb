// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "Poly.hpp"

#include "Range.hpp"
#include "Zip.hpp"
#include <ostream>
#include <iostream>
#include <algorithm>

MATHICGB_NAMESPACE_BEGIN

void Poly::sortTermsDescending() {
  const size_t count = termCount();
  std::vector<size_t> ordered(count);
  for (size_t i = 0; i < count; ++i)
    ordered[i] = i;

  auto cmp = [&](size_t a, size_t b) {
    MATHICGB_ASSERT(a < termCount());
    MATHICGB_ASSERT(b < termCount());
    return monoid().lessThan(mono(b), mono(a));
  };
  std::sort(ordered.begin(), ordered.end(), cmp);

  Poly poly(ring());
  for (size_t i = 0; i < count; ++i)
    poly.append(coef(ordered[i]), mono(ordered[i]));
  *this = std::move(poly);

  MATHICGB_ASSERT(termsAreInDescendingOrder());
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

size_t Poly::getMemoryUse() const {
  size_t total = sizeof(const PolyRing *);
  total += sizeof(coefficient) * mCoefs.capacity();
  total += sizeof(int) * mMonos.capacity();
  return total;
}

bool Poly::termsAreInDescendingOrder() const {
  if (!isZero()) {
    auto prev = leadMono().ptr();
    for (const auto& mono : range(++monoBegin(), monoEnd())) {
      if (monoid().lessThan(*prev, mono))
        return false;
      prev = mono.ptr();
    }
  }
  return true;
}

MATHICGB_NAMESPACE_END
