// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include <iostream>
#include <algorithm>
#include "FreeModuleOrder.hpp"
#include "Poly.hpp"
#include "Ideal.hpp"
#include <mathic.h>
#include "PairTriangle.hpp"

template<class Cmp>
class ConcreteOrder : public FreeModuleOrder {
public:
  ConcreteOrder(Cmp cmp): mCmp(cmp), mComparisons(0), mPreComparisons(0) {}
  virtual ~ConcreteOrder() {}

  virtual int signatureCompare(const_monomial sigA, const_monomial sigB) const {
    ++mComparisons;
    return mCmp.signatureCompare(sigA, sigB);
  }

  virtual int signatureCompare(
    const_monomial sigA,
    const_monomial monoB,
    const_monomial sigB
  ) const {
    ++mComparisons;
    return mCmp.signatureCompare(sigA, monoB, sigB);
  }

  virtual void destructiveSort(std::vector<PreSPair>& pairs) const {
    DestructiveSortComparer cmp(mCmp);
    mCmp.prepareDestructiveSort(pairs);
    std::sort(pairs.begin(), pairs.end(), cmp);
    mPreComparisons += cmp.comparisons();
  }

  virtual void appendBasisElement(const_monomial m) {
    mCmp.appendBasisElement(m);
  }

  virtual std::string description() const {
    return mCmp.description();
  }

  virtual void getStats(size_t& comparisons, size_t& preComparisons) const {
    comparisons = mComparisons;
    preComparisons = mPreComparisons;
  }

private:
  class DestructiveSortComparer {
  public:
    DestructiveSortComparer(const Cmp& cmp): mCmp(cmp), mComparisons(0) {}
    bool operator()(const PreSPair& a, const PreSPair& b) const {
      ++mComparisons;
      return mCmp.lessThanForDestructiveSort(a.signature, b.signature);
    }
    size_t comparisons() const {return mComparisons;}
  private:
    const Cmp& mCmp;
    mutable size_t mComparisons;
  };

  Cmp mCmp;
  mutable size_t mComparisons;
  mutable size_t mPreComparisons;
};

// ** Graded reverse lex.
// Degrees and exponents considered from high index to low index.
//
//   rule: a[i] < b[i] then a > b
//
// also applies to component, which is considered last.
extern int tracingLevel;
class OrderA {
public:
  OrderA(const PolyRing* ring): mRing(ring) {}

  int signatureCompare(const_monomial sigA, const_monomial sigB) const {
    return mRing->monomialCompare(sigA, sigB);
  }

  void prepareDestructiveSort(std::vector<PreSPair>& pairs) const {}

  bool lessThanForDestructiveSort(const_monomial sigA, const_monomial sigB) const {
    return mRing->monomialLT(sigA, sigB);
  }

  int signatureCompare(const_monomial sigA,
    const_monomial monoB, const_monomial sigB) const {
    return mRing->monomialCompare(sigA, monoB, sigB);
  }

  void appendBasisElement(const_monomial) {}

  char const* description() const {return "GrevLex IndexDown";}

private:
  PolyRing const* const mRing;
};

// let d(x) be the first degree of x, so if there are more than one
// vector of weights, only the first one matters for d. By "first"
// degree, read "one with highest index".
//
//   rule 1: if d(ag_i) > (bg_j) then ae_i > be_j
//   rule 2: if i > j then ae_i > be_j
//   rule 2: graded reverse lex as for OrderA, but ignoring the first degree
//
// It was probably intended that all the degree should be considered,
// it just wasn't implemented for more than 1 degree.
class OrderB {
public:
  OrderB(const Ideal *I):
    mRing(I->getPolyRing()),
    topindex(I->getPolyRing()->maxMonomialSize() - 2) {
    for (size_t i = 0; i < I->size(); ++i)
      appendBasisElement(I->getPoly(i)->getLeadMonomial());
  }

  void prepareDestructiveSort(std::vector<PreSPair>& pairs) const {}

  bool lessThanForDestructiveSort(const_monomial a, const_monomial b) const {
    return signatureCompare(a, b) == LT;
  }

  int signatureCompare(const_monomial sig, const_monomial sig2) const {
    int da = - sig[topindex] + deg[*sig];
    int db = -sig2[topindex] + deg[*sig2];
    if (da > db) return GT;
    if (db > da) return LT;
    int cmp = *sig  - *sig2;
    if (cmp > 0) return GT;
    if (cmp < 0) return LT;
    for (size_t i = topindex-1; i >= 1; --i)
      {
        int cmp = sig[i] - sig2[i];
        if (cmp < 0) return GT;
        if (cmp > 0) return LT;
      }
    return EQ;
  }

  int signatureCompare(
    const_monomial sig,
    const_monomial m2,
    const_monomial sig2
  ) const {
    int da = - sig[topindex] + deg[*sig];
    int db = - m2[topindex];
    db += -sig2[topindex] + deg[*sig2];
    if (da > db) return GT;
    if (db > da) return LT;
    int cmp = *sig  - *sig2;
    if (cmp > 0) return GT;
    if (cmp < 0) return LT;
    for (size_t i = topindex-1; i >= 1; --i)
      {
        int cmp = sig[i] - m2[i] - sig2[i];
        if (cmp < 0) return GT;
        if (cmp > 0) return LT;
      }
    return EQ;
  }

  char const* description() const {return "DegreeUp IndexUp GrevLex";}
  void appendBasisElement(const_monomial m) {deg.push_back(-m[topindex]);}

private:
  PolyRing const* const mRing;
  size_t const topindex;

  // array of degrees for each component 0..numgens I - 1
  std::vector<int> deg;
};

// as OrderB, except rule 2 is opposite, so
//   rule 2: if i > j then ae_i < be_j
class OrderC
{
public:
  OrderC(const Ideal* I):
    mRing(I->getPolyRing()),
    topindex(I->getPolyRing()->maxMonomialSize() - 2)
  {
    for (size_t i = 0; i < I->size(); ++i)
      appendBasisElement(I->getPoly(i)->getLeadMonomial());
  }

  void appendBasisElement(const_monomial m) {deg.push_back(-m[topindex]);}

  void prepareDestructiveSort(std::vector<PreSPair>& pairs) const {}

  bool lessThanForDestructiveSort(const_monomial a, const_monomial b) const {
    return signatureCompare(a, b) == LT;
  }

  int signatureCompare(const_monomial sig, const_monomial sig2) const {
    int da = - sig[topindex] + deg[*sig];
    int db = -sig2[topindex] + deg[*sig2];
    if (da > db) return GT;
    if (db > da) return LT;
    int cmp = *sig  - *sig2;
    if (cmp < 0) return GT;
    if (cmp > 0) return LT;
    for (size_t i = topindex-1; i >= 1; --i)
      {
        int cmp = sig[i] - sig2[i];
        if (cmp < 0) return GT;
        if (cmp > 0) return LT;
      }
    return EQ;
  }

  int signatureCompare(
    const_monomial sig,
    const_monomial m2,
    const_monomial sig2
  ) const {
    int da = - sig[topindex] + deg[*sig];
    int db = - m2[topindex];
    db += -sig2[topindex] + deg[*sig2];
    if (da > db) return GT;
    if (db > da) return LT;
    int cmp = *sig  - *sig2;
    if (cmp < 0) return GT;
    if (cmp > 0) return LT;
    for (size_t i = topindex-1; i >= 1; --i)
      {
        int cmp = sig[i] - m2[i] - sig2[i];
        if (cmp < 0) return GT;
        if (cmp > 0) return LT;
      }
    return EQ;
  }

  char const* description() const {return "DegreeUp IndexDown GrevLex";}

private:
  PolyRing const* const mRing;
  const size_t topindex;

  // array of degrees for each component 0..numgens I - 1
  std::vector<int> deg;
};

// Let l(ae_i) be the leading monomial of ag_i.
// Indexes are considered from high to low
// rule 1: higher degree of l(ae_I) wins
// rule 2: reverse lex on l(ae_i)
// rule 3: higher component wins (if mUp is true)
// rule 3': lower components wins (if mUp is false)
class OrderD
{
public:
  OrderD(Ideal const* I, bool up0):
    mRing(I->getPolyRing()),
    mUp(up0),
    topindex(mRing->maxMonomialSize() - 2) {
    for (size_t i = 0; i < I->size(); ++i)
      appendBasisElement(I->getPoly(i)->getLeadMonomial());
  }

  void appendBasisElement(const_monomial m) {monoms.push_back(m);}

  int signatureCompare(const_monomial a, const_monomial b) const {
    const_monomial ma = monoms[*a];
    const_monomial mb = monoms[*b];
    for (size_t i = topindex; i >= 1; i--) {
      int cmp = a[i] - b[i] + ma[i] - mb[i];
      if (cmp != 0)
        return cmp < 0 ? GT : LT;
    }
    if (*a == *b)
      return EQ;
    if (mUp)
      return *a < *b ? LT : GT;
    else
      return *a < *b ? GT : LT;
  }

  void prepareDestructiveSort(std::vector<PreSPair>& pairs) const {
    typedef std::vector<PreSPair>::iterator Iter;
    Iter end = pairs.end();
    for (Iter it = pairs.begin(); it != end; ++it) {
      monomial sig = it->signature;
      const_monomial adjust = monoms[*sig];
      for (size_t i = topindex; i >= 1; --i)
        sig[i] += adjust[i];
      if (mUp)
        *sig = -*sig;
    }
  }

  inline bool lessThanForDestructiveSort(
    const_monomial a,
    const_monomial b
  ) const {
    /* non-unrolled version:
    for (size_t i = topindex; i != static_cast<size_t>(-1); i--) {
      // for i == 0, this will compare the components. We have
      // already adjusted those to take account of mUp.
      int cmp = a[i] - b[i]; // already done: + ma[i] - mb[i];
      if (cmp != 0)
        return cmp > 0;
    }
    return false; // equality */

    size_t i = topindex;
    if ((topindex & 1) == 0) { // if odd number of entries to check
      int cmp = a[topindex] - b[topindex];
      if (cmp != 0)
        return cmp > 0;
      --i;
    }
    for (; i != static_cast<size_t>(-1); i -= 2) {
      int cmp = a[i] - b[i]; // already done: + ma[i] - mb[i];
      if (cmp != 0)
        return cmp > 0;
      int cmp2 = a[i - 1] - b[i - 1];
      if (cmp2 != 0)
        return cmp2 > 0;
    }
    return false; // equality
  }

  int signatureCompare(const_monomial a, const_monomial m2, const_monomial b) const {
    const_monomial ma = monoms[*a];
    const_monomial mb = monoms[*b];
    for (size_t i = topindex; i >= 1; i--)
      {
        int cmp = a[i] - b[i] + ma[i] - mb[i] - m2[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
    int cmp = *a - *b;
    if (mUp)
      {
        if (cmp < 0) return LT;
        if (cmp > 0) return GT;
      }
    else
      {
        if (cmp < 0) return GT;
        if (cmp > 0) return LT;
      }
    return EQ;
  }

  virtual char const* description() const {
    if (mUp)
      return "SchreyerGrevLexUp";
    else
      return "SchreyerGrevLexDown";
  }

private:
  PolyRing const* const mRing;
  bool const mUp; // how to break ties once the monomials are the same
  size_t const topindex; // taken from R

  // array 0..numgens I-1 pointing to lead terms of elements of I
  std::vector<const_monomial> monoms;
};


// Let l(ae_i) be the leading monomial of ag_i.
// Indexes are considered from high to low
// rule 1: higher component wins (if mUp is true)
// rule 1': lower components wins (if mUp is false)
// rule 2: higher degree of l(ae_I) wins
// rule 3: reverse lex on l(ae_i)
class OrderE
{
public:
  OrderE(Ideal const* I, bool up0):
    mRing(I->getPolyRing()),
    mUp(up0),
    topindex(mRing->maxMonomialSize() - 2) {
    for (size_t i = 0; i < I->size(); ++i)
      appendBasisElement(I->getPoly(i)->getLeadMonomial());
  }

  void appendBasisElement(const_monomial m) {monoms.push_back(m);}

  int signatureCompare(const_monomial a, const_monomial b) const {
    const_monomial ma = monoms[*a];
    const_monomial mb = monoms[*b];

    if (*a != *b)
      {
        if (mUp)
          return *a < *b ? LT : GT;
        else
          return *a < *b ? GT : LT;
      }

    for (size_t i = topindex; i >= 1; i--) {
      int cmp = a[i] - b[i] + ma[i] - mb[i];
      if (cmp != 0)
        return cmp < 0 ? GT : LT;
    }
    return EQ;
  }

  void prepareDestructiveSort(std::vector<PreSPair>& pairs) const {
  }

  inline bool lessThanForDestructiveSort(
    const_monomial a,
    const_monomial b
  ) const {
    return signatureCompare(a,b) == LT;
  }

  int signatureCompare(const_monomial a, const_monomial m2, const_monomial b) const {
    const_monomial ma = monoms[*a];
    const_monomial mb = monoms[*b];
    int cmp = *a - *b;
    if (cmp != 0)
      {
        if (mUp)
          {
            if (cmp < 0) return LT;
            if (cmp > 0) return GT;
          }
        else
          {
            if (cmp < 0) return GT;
            if (cmp > 0) return LT;
          }
      }

    for (size_t i = topindex; i >= 1; i--)
      {
        int cmp = a[i] - b[i] + ma[i] - mb[i] - m2[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
    return EQ;
  }

  virtual char const* description() const {
    if (mUp)
      return "IndexUp SchreyerGrevLex";
    else
      return "IndexDown SchreyerGrevLex";
  }

private:
  PolyRing const* const mRing;
  bool const mUp;
  size_t const topindex; // taken from R

  // array 0..numgens I-1 pointing to lead terms of elements of I
  std::vector<const_monomial> monoms;
};

void FreeModuleOrder::displayOrderTypes(std::ostream &o)
{
  o << "FreeModule orders:" << std::endl;
  o << "  1   GrevLex IndexDown" << std::endl;
  o << "  2   DegreeUp IndexUp GrevLex" << std::endl;
  o << "  3   DegreeUp IndexDown GrevLex" << std::endl;
  o << "  4   SchreyerGrevLexUp" << std::endl;
  o << "  5   SchreyerGrevLexDown" << std::endl;
  o << "  6   IndexUp SchreyerGrevLex" << std::endl;
  o << "  7   IndexDown SchreyerGrevLex" << std::endl;
}

FreeModuleOrder *FreeModuleOrder::makeOrder(FreeModuleOrderType type, const Ideal *I)
{
  int i;
  if (type == 0)
    type = 1;  // Set the default

  switch (type) {
  case 1:
    return new ConcreteOrder<OrderA>(OrderA(I->getPolyRing()));
  case 2:
   return new ConcreteOrder<OrderB>(OrderB(I));
  case 3:
   return new ConcreteOrder<OrderC>(OrderC(I));
  case 4:
    return new ConcreteOrder<OrderD>(OrderD(I, true));
  case 5:
    return new ConcreteOrder<OrderD>(OrderD(I, false));
  case 6:
    return new ConcreteOrder<OrderE>(OrderE(I, true));
  case 7:
    return new ConcreteOrder<OrderE>(OrderE(I, false));
  default: break;
  }

  std::cerr << "unknown free module order type" << std::endl;
  std::cerr << "possible orders are: " << std::endl;
  for (i=1; i<=5; i++)
    {
      FreeModuleOrder* F = makeOrder(i,I);
      std::cerr << "  " << i << ": " << F->description() << std::endl;
    }
  exit(1);
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
