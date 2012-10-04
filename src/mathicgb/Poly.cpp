// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include "Poly.hpp"
#include <ostream>
#include <iostream>
#include <algorithm>

// Format for input/output:
//  #terms term1 term2 ...
//  each term: coeff monom
//  each coeff: int
//  each monom: len v1 e1 v2 e2 ... vr er
//   where len = r
// MAJOR ASSUMPTION: the monomials are ordered in descending order!!

void Poly::copy(Poly &result) const
{
  result.R = R;
  result.coeffs.resize(coeffs.size());
  result.monoms.resize(monoms.size());
  std::copy(coeffs.begin(), coeffs.end(), result.coeffs.begin());
  std::copy(monoms.begin(), monoms.end(), result.monoms.begin());
}

void Poly::sortTermsDescending() {
  struct Cmp {
  public:
    Cmp(const Poly& poly): mPoly(poly) {}

    bool operator()(size_t a, size_t b) {
      MATHICGB_ASSERT(a < mPoly.nTerms());
      MATHICGB_ASSERT(b < mPoly.nTerms());
      return mPoly.R->monomialLT(mPoly.monomialAt(b), mPoly.monomialAt(a));
    }

  private:
    const Poly& mPoly;
  };

  const size_t count = nTerms();
  std::vector<size_t> ordered(count);
  for (size_t i = 0; i < count; ++i)
    ordered[i] = i;

  std::sort(ordered.begin(), ordered.end(), Cmp(*this));

  Poly poly(R);
  for (size_t i = 0; i < count; ++i)
    poly.appendTerm(coefficientAt(ordered[i]), monomialAt(ordered[i]));
  *this = std::move(poly);

  MATHICGB_ASSERT(termsAreInDescendingOrder());
}

monomial Poly::monomialAt(size_t index) {
  MATHICGB_ASSERT(index < nTerms());
  return &monoms[index * R->maxMonomialSize()];
}

const_monomial Poly::monomialAt(size_t index) const {
  MATHICGB_ASSERT(index < nTerms());
  return &monoms[index * R->maxMonomialSize()];
}

coefficient& Poly::coefficientAt(size_t index) {
  MATHICGB_ASSERT(index < nTerms());
  return coeffs[index];
}

const coefficient Poly::coefficientAt(size_t index) const {
  MATHICGB_ASSERT(index < nTerms());
  return coeffs[index];
}

void Poly::appendTerm(coefficient a, const_monomial m)
{
  // the monomial will be copied on.
  coeffs.push_back(a);
  size_t len = R->monomialSize(m);
  exponent const * e = m.unsafeGetRepresentation();
  monoms.insert(monoms.end(), e, e + len);
}

void Poly::append(iterator &first, iterator &last)
{
  for ( ; first != last; ++first)
    appendTerm(first.getCoefficient(), first.getMonomial());
}

Poly *Poly::copy() const
{
  Poly *const_this = const_cast<Poly *>(this);
  Poly *result = new Poly(R);
  iterator a = const_this->begin();
  iterator b = const_this->end();
  result->append(a,b);
  return result;
}

void Poly::multByCoefficient(coefficient c)
{
  for (std::vector<coefficient>::iterator i = coeffs.begin(); i != coeffs.end(); i++)
    R->coefficientMultTo(*i, c);
}

void Poly::makeMonic()
{
  if (isZero()) return;
  coefficient c = getLeadCoefficient();
  R->coefficientReciprocalTo(c);
  for (std::vector<coefficient>::iterator i = coeffs.begin(); i != coeffs.end(); i++)
    R->coefficientMultTo(*i, c);
}

bool operator==(const Poly &a, const Poly &b)
{
  const PolyRing* R = a.getRing();
  if (R != b.getRing())
    return false;
  if (a.nTerms() != b.nTerms())
    return false;
  Poly::const_iterator a1 = a.begin();
  Poly::const_iterator b1 = b.begin();
  for ( ; a1 != a.end(); ++a1, ++b1)
    {
      if (a1.getCoefficient() != b1.getCoefficient())
        return false;
      if (!R->monomialEQ(a1.getMonomial(), b1.getMonomial()))
        return false;
    }
  return true;
}

Poly * Poly::add(const PolyRing *R,
                 iterator i,
                 iterator iend,
                 iterator j,
                 iterator jend,
                 size_t &n_compares)
{
  coefficient c;
  n_compares = 0;
  Poly *result = new Poly(R);

  if (i == iend)
    result->append(j, jend);
  else if (j == jend)
    result->append(i, iend);
  else {
    bool done = false;
    while (!done)
    {
      int cmp = R->monomialCompare(i.getMonomial(), j.getMonomial());
      n_compares++;
      switch (cmp) {
      case LT:
        result->appendTerm(j.getCoefficient(), j.getMonomial());
        ++j;
        if (j == jend)
          {
            result->append(i, iend);
            done = true;
          }
        break;
      case GT:
        result->appendTerm(i.getCoefficient(), i.getMonomial());
        ++i;
        if (i == iend)
          {
            result->append(j, jend);
            done = true;
          }
        break;
      case EQ:
        R->coefficientSet(c, i.getCoefficient());
        R->coefficientAddTo(c, j.getCoefficient());
        if (!R->coefficientIsZero(c))
          result->appendTerm(c, i.getMonomial());
        ++j;
        ++i;
        if (j == jend)
          {
            result->append(i, iend);
            done = true;
          }
        else
          {
            if (i == iend)
              {
                result->append(j, jend);
                done = true;
              }
          }
        break;
      }
    }
  }
  return result;
}

const_monomial Poly::backMonomial() const {
  ASSERT(begin() != end());
  return &(monoms.front()) + R->maxMonomialSize() * (nTerms() - 1);
}

void Poly::multByTerm(coefficient a, const_monomial m)
{
  size_t p = 0;
  exponent * n = &monoms[0];
  iterator j = end();

  for (iterator i = begin(); i != j; ++i, ++p, n += R->maxMonomialSize())
    {
      monomial nmon = n;
      R->coefficientMultTo(coeffs[p], a);
      R->monomialMultTo(nmon, m); // changes the monomial pointed to by n.
    }
}
void Poly::multByMonomial(const_monomial m)
{
  size_t p = 0;
  exponent * n = &monoms[0];
  iterator j = end();

  for (iterator i = begin(); i != j; ++i, ++p, n += R->maxMonomialSize())
    {
      monomial nmon = n;
      R->monomialMultTo(nmon, m); // changes the monomial pointed to by n.
    }
}

void Poly::dump() const
{
  std::cout << "coeffs: ";
  for (unsigned i=0; i<coeffs.size(); i++)
    std::cout << " " << coeffs[i];
  std::cout << std::endl;
  std::cout << "monoms: ";
  for (unsigned int i=0; i<monoms.size(); i++)
    std::cout << " " << monoms[i];
  std::cout << std::endl;
}

void Poly::parseDoNotOrder(std::istream& i)
{
  char next = i.peek();
  if (next == '0') {
    i.get();
    return;
  }

  for (;;) {
    bool is_neg = false;
    next = i.peek();
    if (next == '+') {
      i.get();
      next = i.peek();
    } else if (next == '-') {
        is_neg = true;
        i.get();
        next = i.peek();
    } if (isdigit(next) || isalpha(next) || next == '<') {
        int a = 1;
        coefficient b;
        if (isdigit(next))
          {
            i >> a;
            next = i.peek();
          }
        if (is_neg) a = -a;
        size_t firstloc = monoms.size();
        monoms.resize(firstloc + R->maxMonomialSize());
        monomial m = &monoms[firstloc];
        R->coefficientFromInt(b,a);
        coeffs.push_back(b);
        if (isalpha(next) || next == '<')
          R->monomialParse(i, m);
        else
          R->monomialSetIdentity(m); // have to do this to set hash value
        MATHICGB_ASSERT(ring().hashValid(m));
        next = i.peek();
        if (next == '>')
          i.get();
      }
    else
      break;
  }
}

void Poly::parse(std::istream& in) {
  parseDoNotOrder(in);
  sortTermsDescending();
}

void Poly::display(std::ostream &o, bool print_comp) const
{
  long p = R->charac();
  int maxpos = (p+1)/2;
  if (nTerms() == 0)
    {
      o << "0";
      return;
    }
  int nterms = 0;
  for (const_iterator i = begin(); i != end(); ++i)
    {
      nterms++;
      coefficient a = i.getCoefficient();
      bool is_neg = (a > maxpos);
      if (is_neg) a = p-a;
      if (is_neg)
        o << "-";
      else if (nterms > 1)
        o << "+";
      bool print_one = (a == 1);
      if (a != 1) o << a;
      R->monomialDisplay(o, i.getMonomial(), print_comp, print_one);
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
  MATHICGB_ASSERT(R != 0);
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
}

bool Poly::termsAreInDescendingOrder() const {
  if (isZero())
    return true;

  auto stop = end();
  auto it = begin();
  auto previous = it;
  ++it;
  while (it != stop) {
    if (R->monomialCompare(previous.getMonomial(), it.getMonomial()) == LT)
      return false;
    previous = it;
    ++it;
  }
  return true;
}


// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
