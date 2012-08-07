// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include <iostream>
#include "SPairHandler.hpp"
#include "GroebnerBasis.hpp"
#include "MTArray.hpp"
#include "FreeModuleOrder.hpp"
#include "Reducer.hpp"
#include <limits>
#include <stdexcept>

extern int tracingLevel;

SPairHandler::SigPairTriangle::SigPairTriangle(const GroebnerBasis& basis, size_t queueType):
  PairTriangle(basis.order(), basis.ring(), queueType),
  mBasis(basis) {}

bool SPairHandler::SigPairTriangle::calculateOrderBy(
  size_t a,
  size_t b,
  monomial orderBy
) const {
  ASSERT(mBasis.ratioCompare(a, b) != EQ);
  // ensure that ratio(a) > ratio(b)
  if (mBasis.ratioCompare(a, b) == LT)
    std::swap(a, b);
  mBasis.ring().monomialFindSignature(
    mBasis.getLeadMonomial(a),
    mBasis.getLeadMonomial(b),
    mBasis.getSignature(a),
    orderBy);
  return true;
}

SPairHandler::SPairHandler(
  const PolyRing *R0,
  FreeModuleOrder *F0,
  const GroebnerBasis *GB0,
  MonomialTableArray *Hsyz0,
  Reducer* reducer,
  bool postponeKoszuls,
  bool useBaseDivisors,
  bool useSingularCriterionEarly,
  size_t queueType
):
  R(R0),
  F(F0),
  mUseSingularCriterionEarly(useSingularCriterionEarly),
  mUseBaseDivisors(useBaseDivisors),
  mUseHighBaseDivisors(useBaseDivisors),
  Hsyz(Hsyz0),
  GB(GB0),
  mReducer(reducer),
  mPostponeKoszuls(postponeKoszuls),
  mTri(*GB0, queueType) {
}

SPairHandler::~SPairHandler()
{
  ASSERT(mUseBaseDivisors || mUseHighBaseDivisors || mKnownSyzygyTri.empty());
}

void SPairHandler::newSyzygy(const_monomial sig) {
  ASSERT(Hsyz->member(sig));
}

SPairHandler::Stats SPairHandler::getStats() const
{
  F->getStats(mStats.comparisons, mStats.precomparisons);
  return mStats;
}

monomial SPairHandler::popSignature(PairContainer& pairs) {
  ASSERT(!empty());

  monomial sig = R->allocMonomial();
  R->monomialCopy(mTri.topOrderBy(), sig);
  R->setHashOnly(sig); // mTri returns monomials without the hash value set, so we set it here

  pairs.clear();
  do { // pop top as long as the top S-pair has signature sig
    // loop invariant: topGroup is the top of the queue and has signature sig

    { // push the current top S-pair onto pairs
      std::pair<size_t, size_t> p = mTri.topPair();
      size_t greater = p.first;
      size_t smaller = p.second;
      ASSERT(GB->ratioCompare(greater, smaller) != EQ);
      if (GB->ratioCompare(greater, smaller) == LT)
        std::swap(greater, smaller);
      // now, greater in sense of signature in S-pair or, equivalently,
      // in sense of sig/lead ratio.
      ++mStats.spairsFinal;
      pairs.push_back(std::make_pair(greater, smaller));
    }
    mTri.pop();
    ++mStats.duplicateSignatures;
  } while (!mTri.empty() && R->monomialEQ(mTri.topOrderBy(), sig));

  --mStats.duplicateSignatures; // We added one too many in this loop
  return sig;
}

void SPairHandler::newPairs(size_t newGen)
{
  mTri.beginColumn();
  makePreSPairs(newGen);
  mTri.endColumn();

  ASSERT((!mUseBaseDivisors && !mUseHighBaseDivisors) ||
    mKnownSyzygyTri.columnCount() == newGen + 1);
}

void SPairHandler::setupBaseDivisors(
  BaseDivisor& divisor1,
  BaseDivisor& divisor2,
  size_t& highDivisorCmp,
  size_t newGenerator
) {
  BaseDivContainer divisors;
  const size_t MaxBaseDivisors = 2;
  if (mUseBaseDivisors || mUseHighBaseDivisors)
    divisors.reserve(MaxBaseDivisors + 1);

  ASSERT(mUseBaseDivisors || mUseHighBaseDivisors);
  ASSERT(mKnownSyzygyTri.columnCount() == newGenerator);
  mKnownSyzygyTri.addColumn();

  if (mUseHighBaseDivisors) {
    size_t highDivisor = GB->highBaseDivisor(newGenerator);
    if (highDivisor != static_cast<size_t>(-1)) {
      // To use a high divisor, the ratio of the other generator has to be
      // greater than both the ratio of newGenerator and of the high ratio
      // divisor. We can check both at once by letting highDivisorCmp
      // be the one out of newGenerator and highDivisor that has the
      // highest ratio.
      if (GB->ratioCompare(newGenerator, highDivisor) == GT)
        highDivisorCmp = newGenerator;
      else
        highDivisorCmp = highDivisor;
    }
  } else
    highDivisorCmp = static_cast<size_t>(-1);
  if (!mUseBaseDivisors)
    return;

  std::vector<size_t> divs;
  GB->lowBaseDivisors(divs, MaxBaseDivisors, newGenerator);
  ASSERT(divs.size() <= MaxBaseDivisors);

  divisors.resize(divs.size());
  for (size_t i = 0; i < divisors.size(); ++i) {
    BaseDivisor& bd = divisors[i];
    bd.baseDivisor = divs[i];

    // Only use the base divisor technique for generators with ratio
    // less than both N and baseDivisor. baseDivisorCmp is the
    // smallest one of these, so it can be used for this comparison.
    if (GB->ratioCompare(newGenerator, bd.baseDivisor) == LT)
      bd.ratioLessThan = newGenerator;
    else
      bd.ratioLessThan = bd.baseDivisor;

    // Construct a monomial in makeSPair_t2 that can be used
    // to eliminate s-pairs quickly based on the s-pairs already
    // eliminated for baseDivisor.
    const_monomial newSig = GB->getSignature(newGenerator);
    const_monomial newLead = GB->getLeadMonomial(newGenerator);
    const_monomial baseDivSig = GB->getSignature(bd.baseDivisor);
    const_monomial baseDivLead = GB->getLeadMonomial(bd.baseDivisor);
    bd.baseMonomial = R->allocMonomial();
    R->mysteriousSPairMonomialRoutine(newSig, newLead, baseDivSig, baseDivLead, bd.baseMonomial);
  }

  divisor1.baseDivisor = static_cast<size_t>(-1);
  divisor2.baseDivisor = static_cast<size_t>(-1);
  if (divisors.size() >= 1)
    divisor1 = divisors.front();
  if (divisors.size() == 2) {
    divisor2 = divisors.back();
    ASSERT(GB->ratioCompare
      (divisor1.ratioLessThan, divisor2.ratioLessThan) != LT);
  }
}

void SPairHandler::makePreSPairs(size_t newGen)
{
  ASSERT(newGen < GB->size());
  mStats.spairsConstructed += newGen;

  monomial baseDivisorMonomial = 0;

  BaseDivisor divisor1;
  BaseDivisor divisor2;
  divisor1.baseDivisor = static_cast<size_t>(-1);
  divisor2.baseDivisor = static_cast<size_t>(-1);
  size_t highDivisorCmp = static_cast<size_t>(-1);
  if (mUseBaseDivisors || mUseHighBaseDivisors)
    setupBaseDivisors(divisor1, divisor2, highDivisorCmp, newGen);

  monomial hsyz = 0;
  if (!mPostponeKoszuls)
    hsyz = R->allocMonomial();

  const_monomial newSig = GB->getSignature(newGen);
  const_monomial newLead = GB->getLeadMonomial(newGen);
  monomial pairSig = R->allocMonomial();

  if (mUseHighBaseDivisors && divisor1.baseDivisor != static_cast<size_t>(-1))
    ++mStats.hasLowBaseDivisor;
  if (mUseHighBaseDivisors && highDivisorCmp != static_cast<size_t>(-1))
    ++mStats.hasHighBaseDivisor;

  PreSPair result;
  for (size_t oldGen = 0; oldGen < newGen; oldGen++) {
    const_monomial oldSig = GB->getSignature(oldGen);
    const_monomial oldLead = GB->getLeadMonomial(oldGen);

    // Check whether this is a non-regular spair.
    // 'cmp' is used below too.
    const int cmp = GB->ratioCompare(newGen, oldGen);
    if (cmp == EQ) {
      ++mStats.nonregularSPairs;
      continue;
    }

    // check high ratio divisor
    if (mUseHighBaseDivisors &&
      highDivisorCmp != static_cast<size_t>(-1) &&
      GB->ratioCompare(oldGen, highDivisorCmp) == GT &&
      mKnownSyzygyTri.bitUnordered(oldGen, highDivisorCmp)) {
        ASSERT(oldGen != highDivisorCmp); // otherwise ratios should be equal
        mKnownSyzygyTri.setBit(newGen, oldGen, true);
        ++mStats.highBaseDivisorHits;
        // if DEBUG defined, get to the ASSERT below stating
        // that this is really a syzygy
#ifndef DEBUG
        continue;
#endif
    }

    // check low ratio divisors
    if (mUseBaseDivisors &&
      divisor1.baseDivisor != static_cast<size_t>(-1) && 
      GB->ratioCompare(oldGen, divisor1.ratioLessThan) == LT) {
      // if no divisor1, also no divisor 2 and also
      // divisor1 has larger ratio, so skip both checks if divisor1 fails due
      // to the ratio being too small or because there is no divisor1.

      if (
        (divisor1.baseDivisor != oldGen &&  // if divisor1 is a hit
         mKnownSyzygyTri.bitUnordered(divisor1.baseDivisor, oldGen) &&
         R->monomialIsDivisibleBy(divisor1.baseMonomial, oldLead))
        || // or if divisor2 is a hit
        (divisor2.baseDivisor != static_cast<size_t>(-1) && 
         GB->ratioCompare(oldGen, divisor2.ratioLessThan) == LT &&
         divisor2.baseDivisor != oldGen &&
         mKnownSyzygyTri.bitUnordered(divisor2.baseDivisor, oldGen) &&
         R->monomialIsDivisibleBy(divisor2.baseMonomial, oldLead))
      ) {
        mKnownSyzygyTri.setBit(newGen, oldGen, true);
        ++mStats.lowBaseDivisorHits;
        // if DEBUG defined, get to the ASSERT below stating
        // that this really is a syzygy.
#ifndef DEBUG
        continue;
#endif
      }
    }

    if (cmp == GT)
      R->monomialFindSignature(newLead, oldLead, newSig, pairSig);
    else {
      ASSERT(cmp == LT);
      R->monomialFindSignature(oldLead, newLead, oldSig, pairSig);
    }

    size_t result_ignored;
    if (Hsyz->member(pairSig, result_ignored)) {
      ++mStats.syzygyModuleHits;
#ifdef DEBUG
      // Check if actually already elim. by low/high base divisor.
      // Only check in DEBUG mode as otherwise we would have taken an early
      // exit before getting here.
      if ((mUseBaseDivisors || mUseHighBaseDivisors) &&
        mKnownSyzygyTri.bit(newGen, oldGen))
        --mStats.syzygyModuleHits;
#endif
      if (mUseBaseDivisors || mUseHighBaseDivisors)
        mKnownSyzygyTri.setBit(newGen, oldGen, true);
      continue;
    }
    ASSERT((!mUseBaseDivisors && !mUseHighBaseDivisors)
      || !mKnownSyzygyTri.bit(newGen, oldGen));

    if (!mPostponeKoszuls) {
      // add koszul syzygy to Hsyz.
      ASSERT(cmp == GT || cmp == LT);
      if (cmp == GT)
        R->monomialMult(newSig, oldLead, hsyz);
      else
        R->monomialMult(oldSig, newLead, hsyz);
      if (Hsyz->insert(hsyz, 0))
        hsyz = R->allocMonomial();
      if (R->monomialRelativelyPrime(newLead, oldLead))
        {
          ++mStats.earlyRelativelyPrimePairs;
          continue;
        }
    }

    if (mUseSingularCriterionEarly) {
      ASSERT(cmp == GT || cmp == LT);
      size_t const givesSig = (cmp == GT ? newGen : oldGen);    
      if (GB->ratioCompare(GB->minimalLeadInSig(pairSig), givesSig) == GT &&
          !R->monomialRelativelyPrime(newLead, oldLead)) {
        
        ++mStats.earlySingularCriterionPairs;
        continue;
      }
    }

    // construct the PreSPair
    result.signature = pairSig;
    pairSig = R->allocMonomial();
    result.i = static_cast<BigIndex>(oldGen);
    mTri.addPair(result.i, result.signature);
    ++mStats.queuedPairs;
    //pairs.push_back(result);
  }
  R->freeMonomial(pairSig);
  if (mUseBaseDivisors && ! baseDivisorMonomial.isNull())
    R->freeMonomial(baseDivisorMonomial);
  if (!mPostponeKoszuls)
    R->freeMonomial(hsyz);
}

void SPairHandler::setKnownSyzygies(std::vector<std::pair<size_t, size_t> >& pairs) {
  if (!mUseBaseDivisors && !mUseHighBaseDivisors)
    return;
  for (size_t i = 0; i < pairs.size(); ++i)
    setKnownSyzygy(pairs[i].first, pairs[i].second);
}

void SPairHandler::setKnownSyzygy(size_t gen1, size_t gen2) {
  ASSERT(gen1 < GB->size());
  ASSERT(gen2 < GB->size());
  ASSERT(gen1 != gen2);
  if (mUseBaseDivisors || mUseHighBaseDivisors)
    mKnownSyzygyTri.setBitUnordered(gen1, gen2, true);
}

std::string SPairHandler::name() {
  return mTri.name();
}

size_t SPairHandler::getMemoryUse() const
{
  return mTri.getMemoryUse() + getKnownSyzygyBitsMemoryUse();
}

size_t SPairHandler::getKnownSyzygyBitsMemoryUse() const {
  return mKnownSyzygyTri.getMemoryUse();
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
