// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "Poly.hpp"

#include <ostream>
#include <iostream>
#include <algorithm>

MATHICGB_NAMESPACE_BEGIN

// Format for input/output:
//  #terms term1 term2 ...
//  each term: coeff monom
//  each coeff: int
//  each monom: len v1 e1 v2 e2 ... vr er
//   where len = r
// MAJOR ASSUMPTION: the monomials are ordered in descending order!!

Poly& Poly::operator=(Poly&& poly) {
  MATHICGB_ASSERT(&ring() == &poly.ring());
  coeffs = std::move(poly.coeffs);
  monoms = std::move(poly.monoms);
  return *this;
}

void Poly::sortTermsDescending() {
  struct Cmp {
  public:
    Cmp(const Poly& poly): mPoly(poly) {}

    bool operator()(size_t a, size_t b) {
      MATHICGB_ASSERT(a < mPoly.termCount());
      MATHICGB_ASSERT(b < mPoly.termCount());
      return mPoly.ring().monomialLT(mPoly.monomialAt(b), mPoly.monomialAt(a));
    }

  private:
    const Poly& mPoly;
  };

  const size_t count = termCount();
  std::vector<size_t> ordered(count);
  for (size_t i = 0; i < count; ++i)
    ordered[i] = i;

  std::sort(ordered.begin(), ordered.end(), Cmp(*this));

  Poly poly(ring());
  for (size_t i = 0; i < count; ++i)
    poly.appendTerm(coefficientAt(ordered[i]), monomialAt(ordered[i]));
  *this = std::move(poly);

  MATHICGB_ASSERT(termsAreInDescendingOrder());
}

monomial Poly::monomialAt(size_t index) {
  MATHICGB_ASSERT(index < termCount());
  return &monoms[index * ring().maxMonomialSize()];
}

const_monomial Poly::monomialAt(size_t index) const {
  MATHICGB_ASSERT(index < termCount());
  return &monoms[index * ring().maxMonomialSize()];
}

coefficient& Poly::coefficientAt(size_t index) {
  MATHICGB_ASSERT(index < termCount());
  return coeffs[index];
}

const coefficient Poly::coefficientAt(size_t index) const {
  MATHICGB_ASSERT(index < termCount());
  return coeffs[index];
}

void Poly::append(iterator &first, iterator &last)
{
  for ( ; first != last; ++first)
    appendTerm(first.getCoefficient(), first.getMonomial());
}

/*Poly *Poly::copy() const
{
  Poly *const_this = const_cast<Poly *>(this);
  Poly *result = new Poly(*R);
  iterator a = const_this->begin();
  iterator b = const_this->end();
  result->append(a,b);
  return result;
}*/

void Poly::multByCoefficient(coefficient c)
{
  for (std::vector<coefficient>::iterator i = coeffs.begin(); i != coeffs.end(); i++)
    ring().coefficientMultTo(*i, c);
}

bool Poly::isMonic() const {
  return !isZero() && ring().coefficientIsOne(getLeadCoefficient());
}

void Poly::makeMonic() {
  if (isZero())
    return;
  coefficient c = getLeadCoefficient();
  if (ring().coefficientIsOne(c))
    return;
  ring().coefficientReciprocalTo(c);
  for (auto i = coeffs.begin(); i != coeffs.end(); i++)
    ring().coefficientMultTo(*i, c);
  MATHICGB_ASSERT(ring().coefficientIsOne(getLeadCoefficient()));
}

bool operator==(const Poly &a, const Poly &b)
{
  const auto& ring = a.ring();
  if (&a.ring() != &b.ring())
    return false;
  if (a.termCount() != b.termCount())
    return false;
  Poly::const_iterator a1 = a.begin();
  Poly::const_iterator b1 = b.begin();
  for ( ; a1 != a.end(); ++a1, ++b1)
    {
      if (a1.getCoefficient() != b1.getCoefficient())
        return false;
      if (!ring.monomialEQ(a1.getMonomial(), b1.getMonomial()))
        return false;
    }
  return true;
}

const_monomial Poly::backMonomial() const {
  MATHICGB_ASSERT(begin() != end());
  return &(monoms.front()) + ring().maxMonomialSize() * (termCount() - 1);
}

void Poly::multByTerm(coefficient a, const_monomial m)
{
  size_t p = 0;
  exponent * n = &monoms[0];
  iterator j = end();

  for (iterator i = begin(); i != j; ++i, ++p, n += ring().maxMonomialSize())
    {
      monomial nmon = n;
      ring().coefficientMultTo(coeffs[p], a);
      ring().monomialMultTo(nmon, m); // changes the monomial pointed to by n.
    }
}

void Poly::multByMonomial(const_monomial m)
{
  size_t p = 0;
  exponent * n = &monoms[0];
  iterator j = end();

  for (iterator i = begin(); i != j; ++i, ++p, n += ring().maxMonomialSize())
    {
      monomial nmon = n;
      ring().monomialMultTo(nmon, m); // changes the monomial pointed to by n.
    }
}

/*void Poly::dump() const
{
  std::cout << "coeffs: ";
  for (unsigned i=0; i<coeffs.size(); i++)
    std::cout << " " << coeffs[i];
  std::cout << std::endl;
  std::cout << "monoms: ";
  for (unsigned int i=0; i<monoms.size(); i++)
    std::cout << " " << monoms[i];
  std::cout << std::endl;
}*/

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
  const coefficient p = ring().charac();
  const coefficient maxPositive = (p + 1) / 2; // half rounded up
  if (isZero()) {
    out << "0";
    return;
  }
  
  for (const_iterator i = begin(); i != end(); ++i) {
    coefficient coef = i.getCoefficient();
    if (coef > maxPositive) {
      out << "-";
      ring().coefficientNegateTo(coef);
    } else if (i != begin())
      out << '+';
    if (coef != 1)
      out << coef;
    ring().monomialDisplay(out, i.getMonomial(), printComponent, coef == 1);
  }
}

void Poly::display(FILE* file, bool printComponent) const
{
  if (isZero()) {
    fputs("0", file);
    return;
  }

  const auto characteristic = ring().charac();
  const coefficient maxPositiveCoefficient = (characteristic + 1) / 2;
  bool firstTerm = true;
  for (auto it = begin(); it != end(); ++it) {
      coefficient coef = it.getCoefficient();
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
      ring().monomialDisplay(file, it.getMonomial(), printComponent, printOne);
      firstTerm = false;
    }
}

size_t Poly::getMemoryUse() const
{
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
  if (isZero())
    return true;

  auto stop = end();
  auto it = begin();
  auto previous = it;
  ++it;
  while (it != stop) {
    if (ring().monomialCompare(previous.getMonomial(), it.getMonomial()) == LT)
      return false;
    previous = it;
    ++it;
  }
  return true;
}

MATHICGB_NAMESPACE_END
