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

bool PolyRing::hashValid(const_monomial m) const {
  return monomialHashValue(m) == computeHashValue(m);
}

PolyRing::PolyRing(
  coefficient p0,
  int nvars,
  const std::vector<exponent>& weights
):
  mCharac(p0),
  mNumVars(nvars),
  mNumWeights(1),
  mTopIndex(nvars + mNumWeights),
  mHashIndex(nvars + mNumWeights + 1),
  mMaxMonomialSize(nvars + mNumWeights + 2),
  mMaxMonomialByteSize(mMaxMonomialSize * sizeof(exponent)),
  mMonomialPool(mMaxMonomialByteSize),
  mTotalDegreeGradedOnly(false)
#ifdef MATHICGB_USE_MONOID
  , mMonoid(weights)
#endif
#ifdef MATHICGB_USE_FIELD
    , mField(p0)
#endif
{
  MATHICGB_ASSERT(weights.size() == nvars);
  mTotalDegreeGradedOnly = true;
  for (size_t i = 0; i < nvars; ++i)
    if (weights[i] != 1)
      mTotalDegreeGradedOnly = false;
  mWeights.resize(nvars);
  for (size_t i = 0; i < nvars; ++i)
    mWeights[i] = -weights[i];

  resetCoefficientStats();
  srand(0);
  for (size_t i=0; i<mNumVars; i++)
    mHashVals.push_back(static_cast<HashValue>(rand()));
}

PolyRing::PolyRing(coefficient p0,
                   int nvars,
                   int nweights)
  : mCharac(p0),
    mNumVars(nvars),
    mNumWeights(nweights),
    mTopIndex(nvars + nweights),
    mHashIndex(nvars + nweights + 1),
    mMaxMonomialSize(nvars + nweights + 2),
    mMaxMonomialByteSize(mMaxMonomialSize * sizeof(exponent)),
    mMonomialPool(mMaxMonomialByteSize),
    mTotalDegreeGradedOnly(nweights == 1)
#ifdef MATHICGB_USE_MONOID
    , mMonoid(nvars)
#endif
#ifdef MATHICGB_USE_FIELD
    , mField(p0)
#endif
{
  MATHICGB_ASSERT(nweights == 1);

  // set weights to the default value -1
  mWeights.resize(mNumVars * mNumWeights);
  std::fill(mWeights.begin(), mWeights.end(), static_cast<exponent>(-1));

  resetCoefficientStats();
  srand(0);
  for (size_t i=0; i<mNumVars; i++)
    mHashVals.push_back(static_cast<HashValue>(rand()));
}

void PolyRing::displayHashValues() const
{
  std::cerr << "hash values: " << std::endl;
  for (size_t i=0; i<mNumVars; i++)
    std::cerr << "  " << mHashVals[i] << std::endl;
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

bool PolyRing::weightsCorrect(ConstMonomial a1) const
{
  exponent const *a = a1.mValue;
  ++a;
  auto wts = &mWeights[0];
  for (size_t i = 0; i < mNumWeights; ++i) {
    exponent result = 0;
    for (size_t j = 0; j < mNumVars; ++j)
      result += *wts++ * a[j];
    if (a[mNumVars + i] != result)
      return false;
  }
  return true;
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
  for (size_t i = 0; i <= mTopIndex; ++i)
    result[i] = 0;
  result[mHashIndex] = static_cast<HashValue>(0);
}

void PolyRing::monomialEi(size_t i, Monomial &result) const
{
  for (size_t j=mTopIndex; j != static_cast<size_t>(-1); --j)
    result[j] = 0;
  *result  = static_cast<exponent>(i); // todo: handle overflow or change representation
  result[mHashIndex] = static_cast<HashValue>(i); // todo: handle overflow or change representation
}

void PolyRing::monomialMultTo(Monomial &a, ConstMonomial b) const
{
#ifdef MATHICGB_USE_MONOID
  monoid().multiplyInPlace(b, a);
#else
  // a *= b
  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    a[i] += b[i];
#endif
}


void PolyRing::monomialCopy(ConstMonomial a, 
                            Monomial& result) const
{
  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    result[i] = a[i];
}

void PolyRing::monomialQuotientAndMult(ConstMonomial a,
                                       ConstMonomial b,
                                       ConstMonomial c,
                                       Monomial& result) const
{
  // result is set to (a:b) * c.  It is not assumed that b divides a.
  for (size_t i = 0; i <= mNumVars; ++i)
    {
      result[i] = c[i];
      auto cmp = a[i] - b[i];
      if (cmp > 0) result[i] += cmp;
    }
  setWeightsAndHash(result);
}

void PolyRing::monomialFindSignature(ConstMonomial v1,
                                     ConstMonomial v2,
                                     ConstMonomial u1,
                                     Monomial& t1) const
{
  // t1 := (v2:v1) u1
  *t1 = *u1.mValue;
  if (mTotalDegreeGradedOnly) {
    MATHICGB_ASSERT(mNumWeights == 1);
    exponent weight = 0;
    for (size_t i = 1; i <= mNumVars; ++i) {
      MATHICGB_ASSERT(mWeights[i - 1] == -1);
      if (v1[i] < v2[i])
        weight -= t1[i] = u1[i] + v2[i] - v1[i];
      else
        weight -= t1[i] = u1[i];
    }
#ifdef MATHICGB_DEBUG
    setWeightsOnly(t1);
    MATHICGB_ASSERT(t1[mNumVars + 1] == weight);
#endif
    t1[mNumVars + 1] = weight;
  } else {
    for (size_t i = 1; i <= mNumVars; ++i) {
        if (v1[i] < v2[i])
          t1[i] = u1[i] + v2[i] - v1[i];
        else
          t1[i] = u1[i];
    }
    setWeightsOnly(t1);
  }
}

size_t PolyRing::monomialSizeOfSupport(ConstMonomial m) const 
{
  size_t support = 0;
  for (size_t i = 1; i <= mNumVars; ++i)
    if (m[i] != 0)
      ++support;
  return support;
}

void PolyRing::monomialGreatestCommonDivisor(ConstMonomial a, 
                                             ConstMonomial b, 
                                             Monomial& g) const 
{
  *g = 0;
  for (size_t i = 1; i <= mNumVars; ++i)
    g[i] = std::min(a[i], b[i]);
  setWeightsOnly(g);
}

bool PolyRing::monomialIsLeastCommonMultipleNoWeights(
  ConstMonomial a,
  ConstMonomial b,
  ConstMonomial l) const
{
  if (*l != 0)
    return false;
  for (size_t i = 1; i <= mNumVars; ++i)
    if (l[i] != std::max(a[i], b[i]))
      return false;
  return true;
}


void PolyRing::mysteriousSPairMonomialRoutine(ConstMonomial newSig,
                                              ConstMonomial newLead,
                                              ConstMonomial baseDivSig,
                                              ConstMonomial baseDivLead,
                                              Monomial result) const
{
    result[0] = 0; // set component

    for (size_t i = 1; i <= mNumVars; ++i) {
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
  // first initialize result:
  for (size_t j=0; j<mMaxMonomialSize; j++) result[j] = 0;

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
  for (size_t i=0; i<mNumVars; i++) {
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
  for (size_t var = 0; var < mNumVars; ++var) {
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
  for (size_t i = 0; i < mNumVars; ++i) {
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

void PolyRing::printRingFrobbyM2Format(std::ostream& out) const {
  out << "R = QQ[";
  for (size_t i = 0; i < mNumVars; ++i)
    out << (i == 0 ? "" : ", ") << 'x' << (i + 1);
  out << "];\n";
}

PolyRing *PolyRing::read(std::istream &i)
{
  int64 characInt;
  coefficient charac;
  int mNumVars, mNumWeights;

  i >> characInt;
  charac = static_cast<exponent>(characInt);
  i >> mNumVars;
  i >> mNumWeights;
  MATHICGB_ASSERT(mNumWeights == 1);

  std::vector<exponent> weights(mNumWeights);
  int wtlen = mNumVars * mNumWeights;
  weights.resize(wtlen);
  for (int j=0; j < mNumVars * mNumWeights; j++) {
    int64 aInt;
    i >> aInt;
    weights[j] = static_cast<exponent>(aInt);
  }
  return new PolyRing(charac, mNumVars, weights);
}

void PolyRing::write(std::ostream &o) const
{
  o << mCharac << " " << mNumVars << " " << mNumWeights << std::endl;
  auto wts = &mWeights[0];
  for (size_t i = 0; i < mNumWeights; ++i) {
    for (size_t j = 0; j < mNumVars; ++j)
      o << " " << - *wts++;
    o << std::endl;
  }
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
