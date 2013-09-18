// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "Poly.hpp"

#include <ostream>
#include <iostream>
#include <algorithm>

MATHICGB_NAMESPACE_BEGIN

Poly& Poly::operator=(Poly&& poly) {
  MATHICGB_ASSERT(&ring() == &poly.ring());
  coeffs = std::move(poly.coeffs);
  monoms = std::move(poly.monoms);
  return *this;
}

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

auto Poly::mono(size_t index) -> MonoRef {
  MATHICGB_ASSERT(index < termCount());
  return Monoid::toRef(&monoms[index * ring().maxMonomialSize()]);
}

auto Poly::mono(size_t index) const -> ConstMonoRef {
  MATHICGB_ASSERT(index < termCount());
  return Monoid::toRef(&monoms[index * ring().maxMonomialSize()]);
}

coefficient& Poly::coef(size_t index) {
  MATHICGB_ASSERT(index < termCount());
  return coeffs[index];
}

coefficient Poly::coef(size_t index) const {
  MATHICGB_ASSERT(index < termCount());
  return coeffs[index];
}

bool Poly::isMonic() const {
  return !isZero() && ring().coefficientIsOne(leadCoef());
}

void Poly::makeMonic() {
  if (isZero())
    return;
  coefficient c = leadCoef();
  if (ring().coefficientIsOne(c))
    return;
  ring().coefficientReciprocalTo(c);
  for (auto i = coeffs.begin(); i != coeffs.end(); i++)
    ring().coefficientMultTo(*i, c);
  MATHICGB_ASSERT(ring().coefficientIsOne(leadCoef()));
}

bool operator==(const Poly &a, const Poly &b)
{
  if (&a.ring() != &b.ring())
    return false;
  if (a.termCount() != b.termCount())
    return false;
  const auto& ring = a.ring();
  auto a1 = a.begin();
  auto b1 = b.begin();
  for (; a1 != a.end(); ++a1, ++b1) {
    if (a1.coef() != b1.coef())
      return false;
    if (!ring.monoid().equal(a1.mono(), b1.mono()))
      return false;
  }
  return true;
}

const_monomial Poly::backMono() const {
  MATHICGB_ASSERT(begin() != end());
  return &(monoms.front()) + ring().maxMonomialSize() * (termCount() - 1);
}

void Poly::parseDoNotOrder(std::istream& i)
{
  if (i.peek() == '0') {
    i.get();
    return;
  }

  while (true) {
    bool preceededByMinus = false;
    char next = i.peek();
    if (next == '+') {
      i.get();
      next = i.peek();
    } else if (next == '-') {
      preceededByMinus = true;
      i.get();
      next = i.peek();
    }

    if (!isdigit(next) && !isalpha(next) && next != '<')
      break;
    
    { // read coefficient
      int64 bigCoef = 1;
      if (isdigit(next)) {
        i >> bigCoef;
        next = i.peek();
      }
      if (preceededByMinus)
        bigCoef = -bigCoef;
      coeffs.push_back(ring().toCoefficient(bigCoef));
    }

    // read monic monomial
    const size_t firstLocation = monoms.size();
    monoms.resize(firstLocation + ring().maxMonomialSize());
    monomial m = &monoms[firstLocation];
    if (isalpha(next) || next == '<')
      ring().monomialParse(i, m);
    else
      ring().monomialSetIdentity(m); // have to do this to set hash value
    next = i.peek();
    if (next == '>')
      i.get();
  }
}

void Poly::parse(std::istream& in) {
  parseDoNotOrder(in);
  sortTermsDescending();
}

void Poly::display(std::ostream& out, const bool printComponent) const
{
  const auto p = ring().charac();
  const auto maxPositive = (p + 1) / 2; // half rounded up
  if (isZero()) {
    out << "0";
    return;
  }
  
  for (auto i = begin(); i != end(); ++i) {
    auto coef = i.coef();
    if (coef > maxPositive) {
      out << "-";
      ring().coefficientNegateTo(coef);
    } else if (i != begin())
      out << '+';
    if (coef != 1)
      out << coef;
    ring().monomialDisplay(out, Monoid::toOld(i.mono()), printComponent, coef == 1);
  }
}

void Poly::display(FILE* file, bool printComponent) const
{
  if (isZero()) {
    fputs("0", file);
    return;
  }

  const auto characteristic = ring().charac();
  const auto maxPositiveCoefficient = (characteristic + 1) / 2;
  bool firstTerm = true;
  for (auto it = begin(); it != end(); ++it) {
      auto coef = it.coef();
      if (coef > maxPositiveCoefficient) {
        coef = characteristic - coef;
        fputc('-', file);
      } else if (!firstTerm)
        fputc('+', file);
      bool printOne = true;
      if (coef != 1) {
        printOne = false;
        fprintf(file, "%li", (long)coef);
      }
      ring().monomialDisplay(file, Monoid::toOld(it.mono()), printComponent, printOne);
      firstTerm = false;
    }
}

size_t Poly::getMemoryUse() const {
  size_t total = sizeof(const PolyRing *);
  total += sizeof(coefficient) * coeffs.capacity();
  total += sizeof(int) * monoms.capacity();
  return total;
}

void Poly::setToZero() {
  coeffs.clear();
  monoms.clear();
}

void Poly::see(bool print_comp) const
{
  display(std::cout, print_comp);
  std::cout << std::endl;
}

std::ostream& operator<<(std::ostream& out, const Poly& p) {
  p.see(false);
  return out;
}

void Poly::reserve(size_t spaceForThisManyTerms) {
  monoms.reserve(spaceForThisManyTerms * ring().maxMonomialSize());
}

bool Poly::termsAreInDescendingOrder() const {
  if (!isZero()) {
    auto prev = leadMono().ptr();
    MonoRange range = {++monoBegin(), monoEnd()};
    for (const auto& mono : range) {
      if (monoid().lessThan(*prev, mono))
        return false;
      prev = mono.ptr();
    }
  }
  return true;
}

MATHICGB_NAMESPACE_END
