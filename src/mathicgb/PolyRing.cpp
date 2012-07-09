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

PolyRing::PolyRing(long p0,
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
    mTotalDegreeGradedOnly(false)
{
  resetCoefficientStats();
  //  rand();
  for (size_t i=0; i<mNumVars; i++)
    mHashVals.push_back(rand());
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

inline bool PolyRing::weightsCorrect(ConstMonomial a1) const
{
  exponent const *a = a1.mValue;
  ++a;
  const int *wts = &mWeights[0];
  for (int i=0; i < mNumWeights; ++i) {
    int result = 0;
    for (size_t j=0; j < mNumVars; ++j)
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
      int cmp = sig[i] - m2[i] - sig2[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
  return EQ;
}

void PolyRing::monomialSetIdentity(Monomial& result) const
{
  for (size_t i = 0; i <= mTopIndex; ++i)
    result[i] = 0;
  result[mHashIndex] = 0;
}

void PolyRing::monomialEi(size_t i, Monomial &result) const
{
  for (size_t j=mTopIndex; j != static_cast<size_t>(-1); --j) result[j] = 0;
  *result  = static_cast<int>(i); // todo: handle overflow or change representation
  result[mHashIndex] = static_cast<int>(i); // todo: handle overflow or change representation
}

void PolyRing::monomialMultTo(Monomial &a, ConstMonomial b) const
{
  // a *= b
  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    a[i] += b[i];
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
      int cmp = a[i] - b[i];
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
    ASSERT(mNumWeights == 1);
    exponent weight = 0;
    for (size_t i = 1; i <= mNumVars; ++i) {
      ASSERT(mWeights[i - 1] == -1);
      if (v1[i] < v2[i])
        weight -= t1[i] = u1[i] + v2[i] - v1[i];
      else
        weight -= t1[i] = u1[i];
    }
#ifdef DEBUG
    setWeightsOnly(t1);
    ASSERT(t1[mNumVars + 1] == weight);
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

bool PolyRing::monomialIsLeastCommonMultiple(
  ConstMonomial a,
  ConstMonomial b,
  ConstMonomial l) const
{
  return monomialIsLeastCommonMultipleNoWeights(a, b, l) && weightsCorrect(l);
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

  int v, e, x;
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
            {
              i >> e;
            }
          result[v+1] = e;
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
  for (size_t i=0; i<mNumVars; i++)
    if (a[i+1] != 0)
      {
        printed_any = true;
        if (i <= 25)
          o << static_cast<unsigned char>('a' + i);
        else
          o << static_cast<unsigned char>('A' + (i - 26));
        if (a[i+1] != 1)
          o << a[i+1];
      }
  if (print_comp)
    o << '<' << *a.mValue << '>';
  else if (!printed_any && print_one)
    o << "1";
}

void PolyRing::printMonomialFrobbyM2Format(std::ostream& out, ConstMonomial m) const {
  out << "  ";
  bool isOne = true;
  for (size_t i = 0; i < mNumVars; ++i) {
    int e = m[i + 1];
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



///////////////////////////////////////
// (Old) Monomial Routines ////////////
///////////////////////////////////////

#ifndef NEWMONOMIALS
bool PolyRing::monomialIsLeastCommonMultiple(
  const_monomial a,
  const_monomial b,
  const_monomial l) const
{
  return monomialIsLeastCommonMultipleNoWeights(a, b, l) && weightsCorrect(l);
}

inline bool PolyRing::weightsCorrect(const_monomial a) const
{
  ++a;
  const int *wts = &mWeights[0];
  for (int i=0; i < mNumWeights; ++i) {
    int result = 0;
    for (size_t j=0; j < mNumVars; ++j)
      result += *wts++ * a[j];
    if (a[mNumVars + i] != result)
      return false;
  }
  return true;
}

bool PolyRing::monomialIsLeastCommonMultipleNoWeights(
  const_monomial a,
  const_monomial b,
  const_monomial l) const
{
  if (*l != 0)
    return false;
  for (size_t i = 1; i <= mNumVars; ++i)
    if (l[i] != std::max(a[i], b[i]))
      return false;
  return true;
}

int PolyRing::monomialCompare(const_monomial sig, const_monomial m2, const_monomial sig2) const
  // returns LT, EQ, or GT, depending on sig ? (m2 * sig2).
{
  for (size_t i = mTopIndex; i != static_cast<size_t>(-1); --i)
    {
      int cmp = sig[i] - m2[i] - sig2[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
  return EQ;
}

bool PolyRing::monomialEQ(const_monomial a, const_monomial b) const
{
  for (size_t i = 0; i <= mNumVars; ++i)
    if (a[i] != b[i]) return false;
  return true;
}

void PolyRing::monomialSetIdentity(monomial& result) const {
  for (size_t i = 0; i <= mTopIndex; ++i)
    result[i] = 0;
  result[mHashIndex] = 0;
}

void PolyRing::monomialEi(size_t i, monomial &result) const
{
  for (size_t j=mTopIndex; j != static_cast<size_t>(-1); --j) result[j] = 0;
  *result = static_cast<int>(i); // todo: handle overflow or change representation
  result[mHashIndex] = static_cast<int>(i); // todo: handle overflow or change representation
}

void PolyRing::monomialMult(const_monomial a, const_monomial b, monomial &result) const
{
  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    result[i] = a[i] + b[i];
}
void PolyRing::monomialMultTo(monomial a, const_monomial b) const
{
  // a *= b
  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    a[i] += b[i];
}


void PolyRing::monomialCopy(const_monomial a, monomial &result) const
{
  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    result[i] = a[i];
}
void PolyRing::monomialQuotientAndMult(const_monomial a,
                                       const_monomial b,
                                       const_monomial c,
                                       monomial result) const
{
  // result is set to (a:b) * c.  It is not assumed that b divides a.
  for (size_t i = 0; i <= mNumVars; ++i)
    {
      result[i] = c[i];
      int cmp = a[i] - b[i];
      if (cmp > 0) result[i] += cmp;
    }
  setWeightsAndHash(result);
}

void PolyRing::setWeightsAndHash(monomial a) const
{
  setWeightsOnly(a);
  int hash = *a;
  a++;
  for (size_t i = 0; i < mNumVars; ++i)
    hash += a[i] * mHashVals[i];
  a[mHashIndex - 1] = hash;
}

void PolyRing::monomialFindSignatures(const_monomial v1,
                                      const_monomial v2,
                                      const_monomial u1,
                                      const_monomial u2,
                                      monomial t1,
                                      monomial t2) const
{
  // t1 := (v2:v1) u1
  // t2 := (v1:v2) u2
  // set components:
  *t1 = *u1;
  *t2 = *u2;

  for (size_t i = 1; i <= mNumVars; ++i)
    {
      exponent a = v1[i];
      exponent b = v2[i];
      exponent min = std::min(a, b);
      t1[i] = b - min + u1[i];
      t2[i] = a - min + u2[i];
    }
  setWeightsOnly(t1);
  setWeightsOnly(t2);
}

void PolyRing::monomialFindSignature(const_monomial v1,
                                     const_monomial v2,
                                     const_monomial u1,
                                     monomial t1) const
{
  // t1 := (v2:v1) u1
  *t1 = *u1;
  if (mTotalDegreeGradedOnly) {
    ASSERT(mNumWeights == 1);
    exponent weight = 0;
    for (size_t i = 1; i <= mNumVars; ++i) {
      ASSERT(mWeights[i - 1] == -1);
      if (v1[i] < v2[i])
        weight -= t1[i] = u1[i] + v2[i] - v1[i];
      else
        weight -= t1[i] = u1[i];
    }
#ifdef DEBUG
    setWeightsOnly(t1);
    ASSERT(t1[mNumVars + 1] == weight);
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

size_t PolyRing::monomialSizeOfSupport(const_monomial m) const {
  size_t support = 0;
  for (size_t i = 1; i <= mNumVars; ++i)
    if (m[i] != 0)
      ++support;
  return support;
}

void PolyRing::monomialGreatestCommonDivisor(const_monomial a, const_monomial b, monomial g) const {
  *g = 0;
  for (size_t i = 1; i <= mNumVars; ++i)
    g[i] = std::min(a[i], b[i]);
  setWeightsOnly(g);
}

void PolyRing::monomialRead(std::istream &I, monomial &result) const
{
  // format is as stated in PolyRing.h
  int len;
  I >> *result; // component
  I >> len; // nterms of monomial
  for (size_t i=1; i<=mNumVars; i++)
    result[i] = 0;
  for (int i=len; i>0; --i)
    {
      int v,e;
      I >> v;
      I >> e;
      result[v+1] = e;
    }
  setWeightsAndHash(result);
}
void PolyRing::monomialWrite(std::ostream &o, const_monomial monom) const
{
  int len = 0;
  for (size_t i=1; i<=mNumVars; i++)
    if (monom[i] != 0) len++;
  o << *monom << " " << len; // the component
  monom++;
  for (size_t i = mNumVars - 1; i != static_cast<size_t>(-1); --i)
    if (monom[i] != 0)
      o << " " << i << " " << monom[i];
}

// The following two read/write monomials in the form:
//  letter number letter number ... letter number < number >
// allowed letters: a-zA-Z
// also allowed instead of letter: [ number ], [0] refers to first var, etc.
// What about errors??
void PolyRing::monomialParse(std::istream &i, monomial &result) const
{
  // first initialize result:
  for (size_t j=0; j<mMaxMonomialSize; j++) result[j] = 0;

  int v, e, x;
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
            {
              i >> e;
            }
          result[v+1] = e;
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

void PolyRing::monomialDisplay(std::ostream &o, const_monomial a, bool print_comp, bool print_one) const
{
  // We display a monomial in the form that can be used with monomialParse
  // print_one: only is consulted in print_comp is false.

  bool printed_any = false;
  for (size_t i=0; i<mNumVars; i++)
    if (a[i+1] > 0)
      {
        printed_any = true;
        if (i <= 25)
          o << static_cast<unsigned char>('a' + i);
        else
          o << static_cast<unsigned char>('A' + (i - 26));
        if (a[i+1] >= 2)
          o << a[i+1];
      }
  if (print_comp)
    o << '<' << *a << '>';
  else if (!printed_any && print_one)
    o << "1";
}

void PolyRing::printMonomialFrobbyM2Format(std::ostream& out, const_monomial m) const {
  out << "  ";
  bool isOne = true;
  for (size_t i = 0; i < mNumVars; ++i) {
    int e = m[i + 1];
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
#endif

PolyRing *PolyRing::read(std::istream &i)
{
  long charac;
  int mNumVars, mNumWeights;

  i >> charac;
  i >> mNumVars;
  i >> mNumWeights;
  PolyRing *R = new PolyRing(charac, mNumVars, mNumWeights);
  int wtlen = mNumVars * mNumWeights;
  R->mWeights.resize(wtlen);
  R->mTotalDegreeGradedOnly = (mNumWeights == 1);
  for (int j=0; j <mNumVars * mNumWeights; j++)
    {
      int a;
      i >> a;
      R->mWeights[j] = -a;
      if (R->mWeights[j] != -1)
        R->mTotalDegreeGradedOnly = false;
    }
  return R;
}

void PolyRing::write(std::ostream &o) const
{
  o << mCharac << " " << mNumVars << " " << mNumWeights << std::endl;
  const int *wts = &mWeights[0];
  for (int i=0; i<mNumWeights; i++)
    {
      for (size_t j=0; j<mNumVars; j++)
        o << " " << - *wts++;
      o << std::endl;
    }
}

void PolyRing::coefficientFromInt(coefficient &result, int a) const
{
  result = a % mCharac;
  if (result < 0) result += mCharac;
}

void PolyRing::coefficientAddOneTo(coefficient &result) const
{
  result++;
  if (result == mCharac) result = 0;
}
void PolyRing::coefficientNegateTo(coefficient &result) const
 // result = -result
{
  result = mCharac - result;
}
void PolyRing::coefficientAddTo(coefficient &result, coefficient a, coefficient b) const
// result += a*b
{
  mStats.n_addmult++;
  long c = a * b + result;
  result = c % mCharac;
}
void PolyRing::coefficientAddTo(coefficient &result, coefficient a) const
 // result += a
{
  mStats.n_add++;
  result += a;
  if (result >= mCharac) result -= mCharac;
}
void PolyRing::coefficientMultTo(coefficient &result, coefficient a) const
  // result *= a
{
  mStats.n_mult++;
  long b = result * a;
  result = b % mCharac;
}
void PolyRing::coefficientMult(coefficient a, coefficient b, coefficient &result) const
{
  mStats.n_mult++;
  long c = b * a;
  result = c % mCharac;
}

void gcd_extended(long a, long b, long &u, long &v, long &g)
{
  long q ;
  long u1, v1, g1;
  long utemp, vtemp, gtemp;

  g1 = b;     u1 = 0;         v1 = 1;
   g = a;      u = 1;          v = 0;
  while (g1 != 0)
    {
      q = g / g1 ;
      gtemp= g - q * g1;
      utemp= u - q * u1;
      vtemp= v - q * v1;
       g = g1;             u = u1;                 v = v1 ;
      g1 = gtemp;         u1 = utemp;             v1 = vtemp;
    }
}

void PolyRing::coefficientReciprocalTo(coefficient &result) const
{
  mStats.n_recip++;
  long u,v,g;
  gcd_extended(result, mCharac, u, v, g);
  if (u < 0) u += mCharac;
  result = u;
}

void PolyRing::coefficientDivide(coefficient a, coefficient b, coefficient &result) const
 // result /= a
{
  mStats.n_divide++;
  long u,v,g;
  gcd_extended(b, mCharac, u, v, g);
  result = a * u;
  result %= mCharac;
  if (result < 0) result += mCharac;
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
