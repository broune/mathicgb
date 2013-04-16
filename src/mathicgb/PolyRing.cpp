// Copyright 2011 Michael E. Stillman
#include "stdinc.h"
#include "PolyRing.hpp"

#include <algorithm>
#include <ostream>
#include <istream>
#include <iostream>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <limits>

PolyRing::PolyRing(const Field& field, const Monoid& monoid):
  mField(field),
  mMonoid(monoid),
  mNumWeights(monoid.gradingCount()),
  mTopIndex(monoid.varCount() + mNumWeights),
  mMonomialPool(maxMonomialByteSize())
{
}

PolyRing::PolyRing(
  coefficient p0,
  int nvars,
  const std::vector<exponent>& weights
):
  mField(p0),
  mMonoid(nvars, weights),
  mNumWeights(nvars == 0 ? 0 : weights.size() / nvars),
  mTopIndex(nvars + mNumWeights),
  mMonomialPool(maxMonomialByteSize())
{
}

PolyRing::PolyRing(coefficient p0,
                   int nvars,
                   int nweights):
  mField(p0),
  mMonoid(nvars),
  mNumWeights(nweights),
  mTopIndex(nvars + nweights),
  mMonomialPool(maxMonomialByteSize())
{
}

void PolyRing::resetCoefficientStats() const
{
  mStats.n_addmult = 0;
  mStats.n_add = 0;
  mStats.n_mult = 0;
  mStats.n_recip = 0;
  mStats.n_divide = 0;
}

///////////////////////////////////////
// (New) Monomial Routines ////////////
///////////////////////////////////////

void PolyRing::setWeightsAndHash(Monomial& a1) const
{
  setWeightsOnly(a1);
  setHashOnly(a1);
}

int PolyRing::monomialCompare(ConstMonomial sig, ConstMonomial m2, ConstMonomial sig2) const
  // returns LT, EQ, or GT, depending on sig ? (m2 * sig2).
{
  for (size_t i = mTopIndex; i != static_cast<size_t>(-1); --i)
    {
      auto cmp = sig[i] - m2[i] - sig2[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
  return EQ;
}

void PolyRing::monomialSetIdentity(Monomial& result) const
{
  monoid().setIdentity(result);
}

void PolyRing::monomialEi(size_t i, Monomial &result) const
{
  monoid().setIdentity(result);
  monoid().setComponent(i, result);
}

void PolyRing::monomialMultTo(Monomial &a, ConstMonomial b) const
{
  monoid().multiplyInPlace(b, a);
}


void PolyRing::monomialCopy(ConstMonomial a, 
                            Monomial& result) const
{
  monoid().copy(a, result);
}

void PolyRing::monomialQuotientAndMult(ConstMonomial a,
                                       ConstMonomial b,
                                       ConstMonomial c,
                                       Monomial& result) const
{
  monoid().colonMultiply(b, a, c, result);
}

void PolyRing::monomialFindSignature(ConstMonomial v1,
                                     ConstMonomial v2,
                                     ConstMonomial u1,
                                     Monomial& t1) const
{
  monoid().colonMultiply(v1, v2, u1, t1);
}

size_t PolyRing::monomialSizeOfSupport(ConstMonomial m) const 
{
  return monoid().sizeOfSupport(m);
}

void PolyRing::mysteriousSPairMonomialRoutine(ConstMonomial newSig,
                                              ConstMonomial newLead,
                                              ConstMonomial baseDivSig,
                                              ConstMonomial baseDivLead,
                                              Monomial result) const
{
    result[0] = 0; // set component
    for (size_t i = 1; i <= varCount(); ++i) {
      exponent x = newSig[i] - baseDivSig[i] + baseDivLead[i];
      if (newLead[i] <= x)
        x = std::numeric_limits<exponent>::max();
      result[i] = std::max(baseDivLead[i], x);
    }
}

// The following two read/write monomials in the form:
//  letter number letter number ... letter number < number >
// allowed letters: a-zA-Z
// also allowed instead of letter: [ number ], [0] refers to first var, etc.
// What about errors??
void PolyRing::monomialParse(std::istream &i, 
                             Monomial& result) const
{
//  monoid().parseM2(i, result);
  // first initialize result:
  for (size_t j=0; j< maxMonomialSize(); j++) result[j] = 0;

  uint64 e;
  int v, x;
  // now look at the next char
  for (;;)
    {
      char next = i.peek();
      if (isalpha(next))
        {
          // determine variable
          v = next - 'a';
          if (v > 25 || v < 0)
            v = next - 'A' + 26;
          if (v < 0 || v > 51)
            std::cerr << "variable out of range" << std::endl;
          i.get(); // move past the letter
          next = i.peek();
          e = 1;
          if (isdigit(next))
            i >> e;
          result[v+1] = static_cast<exponent>(e);
        }
      else if (next == '<')
        {
          // get component
          i.get();
          i >> x;
          *result = x; // place component
          if (i.peek() == '>') i.get();
        }
      else
        {
          break; // We assume that we are at the end of the monomial
        }
    }
  setWeightsAndHash(result);
}

void PolyRing::monomialDisplay(std::ostream &o, 
                               ConstMonomial a, 
                               bool print_comp, 
                               bool print_one) const
{
  // We display a monomial in the form that can be used with monomialParse
  // print_one: only is consulted in print_comp is false.

  bool printed_any = false;
  for (size_t i=0; i< varCount(); i++) {
    if (a[i+1] != 0) {
      printed_any = true;
      if (i <= 25)
        o << static_cast<unsigned char>('a' + i);
      else
        o << static_cast<unsigned char>('A' + (i - 26));
      if (a[i+1] != 1)
        o << static_cast<int64>(a[i+1]);
    }
  }
  if (print_comp)
    o << '<' << static_cast<int64>(*a.mValue) << '>';
  else if (!printed_any && print_one)
    o << "1";
}

void PolyRing::monomialDisplay(FILE* file, 
                               ConstMonomial mono, 
                               bool printComponent, 
                               bool printOne) const
{
  const unsigned int letterCount = 26;
  MATHICGB_ASSERT(getNumVars() <= 2 * letterCount);
  bool printedAny = false;
  for (size_t var = 0; var < varCount(); ++var) {
    exponent e = monomialExponent(mono, var);
    if (e == 0)
      continue;
    printedAny = true;
    char varChar;
    if (var < letterCount)
      varChar = static_cast<char>('a' + var);
    else
      varChar = static_cast<char>('A' + (var - letterCount));
    fputc(varChar, file);
    if (e != 1)
      fprintf(file, "%i", e);
  }
  if (printComponent)
    fprintf(file, "<%i>", mono.component());
  else if (!printedAny && printOne)
    fputc('1', file);
}

void PolyRing::printMonomialFrobbyM2Format(std::ostream& out, ConstMonomial m) const {
  out << "  ";
  bool isOne = true;
  for (size_t i = 0; i < varCount(); ++i) {
    const auto e = m[i + 1];
    if (e == 0)
      continue;
    if (!isOne)
      out << '*';
    else
      isOne = false;
    out << 'x' << (i + 1) << '^' << e;
  }
  if (isOne)
    out << '1';
}

PolyRing *PolyRing::read(std::istream &i)
{
  int64 characInt;
  coefficient charac;
  i >> characInt;
  charac = static_cast<exponent>(characInt);
  return new PolyRing(charac, Monoid::readMonoid(i));
}

void PolyRing::write(std::ostream &o) const
{
  o << charac() << ' ';
  monoid().printMonoid(o);
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
