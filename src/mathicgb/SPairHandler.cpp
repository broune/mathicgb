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
  mTrackEssentialPair(false),
  mEssentialFirst(static_cast<size_t>(-1)),
  mEssentialSecond(static_cast<size_t>(-1)),
  mEssentialSig(R0->allocMonomial()),
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
  R->freeMonomial(mEssentialSig);
}

bool SPairHandler::hasEssentialPair() const { 
  ASSERT(mTrackEssentialPair);
  bool value = (mEssentialFirst < GB->size());
  ASSERT(!value || isEssential(mEssentialFirst, mEssentialSecond));
  return value;
}

void SPairHandler::newSyzygy(const_monomial sig) {
  ASSERT(Hsyz->member(sig));
}

bool SPairHandler::isEssential(size_t a, size_t b) const {
  ASSERT(mTrackEssentialPair);
  bool minA = GB->basis().leadMinimal(a);
  bool minB = GB->basis().leadMinimal(b);
  if (minA && minB)
    return true;
  if (!minA && !minB)
    return false;
  if (!minA)
    std::swap(a, b);
  ASSERT(GB->basis().leadMinimal(a));
  ASSERT(!GB->basis().leadMinimal(b));
  if (mDidReducingSPair[b])
    return false;
  return R->monomialIsDivisibleBy
    (GB->getLeadMonomial(b), GB->getLeadMonomial(a));
}

void SPairHandler::setTrackEssentialPair(bool value) {
  if (value == mTrackEssentialPair)
    return;
  mTrackEssentialPair = value;
  if (mTrackEssentialPair) {
    ASSERT(mEssentialPoly.get() == 0);
    while (mDidReducingSPair.size() < GB->size())
      mDidReducingSPair.push_back(false);
    mEssentialFirst = 0;
    mEssentialSecond = static_cast<size_t>(-1);
    if (GB->size() > 0)
      nextEssentialPair();
  } else {
    mEssentialPoly.reset(0);
    mEssentialFirst = static_cast<size_t>(-1);
    mEssentialSecond = static_cast<size_t>(-1);
  }
}

void SPairHandler::nextEssentialPair() {
  ASSERT(mTrackEssentialPair);
  ASSERT(mEssentialFirst < GB->size());
  ASSERT(mEssentialSecond == static_cast<size_t>(-1) ||
    mEssentialSecond < mEssentialFirst);

  if (mEssentialPoly.get() != 0)
    mEssentialPoly.reset(0);

  // If the queue is empty, we don't know what the next signature will be,
  // so we don't know what has been shown to reduce to zero, which means
  // we can't proceed. The queue can be empty without the computation
  // being done if the last S-pair does not reduce to zero, and this
  // method gets called between the last S-pair being popped and the
  // new basis element being inserted.
  ASSERT(!mTri.empty());
  const_monomial currentSig = mTri.topOrderBy();
  const size_t basisElementCount = GB->size();

  size_t& first = mEssentialFirst;
  size_t& second = mEssentialSecond;
  bool firstIteration = true;
  while (true) {
    if (!mDidReducingSPair[first] &&
      !firstIteration &&
      !GB->basis().leadMinimal(first)) {
      ASSERT(R->monomialIsDivisibleBy
        (GB->getLeadMonomial(first), GB->getLeadMonomial(second)));
      mDidReducingSPair[first] = true;
      second = first; // move on to next first
    } else {
nextIteration:
      firstIteration = false;
      ++second;
    }
    // ** move on to next S-pair
    // Setting mEssentialSecond to -1 instructs this method to search the
    // S-pairs with the given mEssentialFirst starting at the first one.
    if (first == second) {
      ++first;
      if (first == basisElementCount) {
        second = static_cast<size_t>(-1); // start at zero next time
        break; // no more S-pairs to check for now
      }
      second = 0;
    }

    // Let
    //   A = the set of basis elements with minimal lead term
    //   B = the set of basis elements with non-minimal lead term
    //   G = the union of A and B
    //   P' = the set of pairs (a,b) in AxB such that lead(a)|lead(b)
    //   P = The union of AxA and P'
    // We say that an S-pair is essential if it is a member of P.

    // Skip pair if not a member of P since we only need the S-polynomials of
    // members of P to reduce to zero. We can skip all the other S-pairs in GxG
    // due to Buchberger's second criterion (lcm criterion).

    if (!GB->basis().leadMinimal(first)) {
      if (mDidReducingSPair[first]) {
        second = first - 1;
        continue;
      }
      for (; second < first; ++second) {
        if (GB->basis().leadMinimal(second) &&
          R->monomialIsDivisibleBy
          (GB->getLeadMonomial(first), GB->getLeadMonomial(second))) {
          ASSERT(isEssential(first, second));
          break;
        }
        ASSERT(!isEssential(first, second));
      }
    } else {
      for (; second < first; ++second) {
        if (GB->basis().leadMinimal(second))
          break;
        if (!mDidReducingSPair[second] &&
          R->monomialIsDivisibleBy
            (GB->getLeadMonomial(second), GB->getLeadMonomial(first))) {
          ASSERT(isEssential(first, second));
          break;
        }
        ASSERT(!isEssential(first, second));
      }
    }

    if (second == first) {
      --second; // as it is incremented at top of loop
      goto nextIteration; // todo: we really need to improve the code structure here
    }
    ASSERT(isEssential(first, second));

    ASSERT(mTrackEssentialPair);
    bool minFirst = GB->basis().leadMinimal(first);
    bool minSecond = GB->basis().leadMinimal(second);
    if (minFirst) {
      if (!minSecond) {
        if (!R->monomialIsDivisibleBy
          (GB->getLeadMonomial(second), GB->getLeadMonomial(first))) {
          ASSERT(!isEssential(first, second));
          continue;
        }
      }
    } else {
      if (!minSecond) {
        ASSERT(!isEssential(first, second));
        continue;
      }
      if (!R->monomialIsDivisibleBy
        (GB->getLeadMonomial(first), GB->getLeadMonomial(second))) {
        ASSERT(!isEssential(first, second));
        continue;
      }
    }
    ASSERT(isEssential(first, second));
    //if (!isEssential(first, second))
    //  continue;

    // Use buchberger's first criterion: skip S-pairs whose leading terms
    // are relatively prime.
    if (R->monomialRelativelyPrime(
      GB->getLeadMonomial(first),
      GB->getLeadMonomial(second)))
      continue;

    // Use Buchberger's second criterion (lcm criterion).
    if (GB->basis().buchbergerLcmCriterion(first, second))
      continue;

    // Skip pair if signature less than or equal to current signature, as
    // then the S-polynomial has been reduced to zero. Note that it is not
    // relevant if the signature of the S-pair is an element of the initial
    // syzygy submodule, as that only implies that the S-pair will reduce to
    // once we get to its signature. It does not imply that the S-pair reduces
    // to zero right now.
    //
    // Usually in the Signature Buchberger Algorithm, it is not OK for
    // the signature of the two halves of the S-pair to be equal.
    // So those S-pairs are usually discarded. We can't do that here.
    size_t greater = first;
    size_t smaller = second;
    if (GB->ratioCompare(greater, smaller) == LT)
      std::swap(greater, smaller);
    R->monomialFindSignature(
      GB->getLeadMonomial(greater),
      GB->getLeadMonomial(smaller),
      GB->getSignature(greater),
      mEssentialSig); // Set mEssentialSig to signature of S-pair

    if (F->signatureCompare(mEssentialSig, currentSig) == LT)
      continue;

    mEssentialPoly = mReducer->classicReduceSPoly
      (GB->poly(first), GB->poly(second), GB->basis());
    ASSERT(mEssentialPoly.get() != 0);
    if (mEssentialPoly->isZero()) {
      mEssentialPoly.reset(0);
      continue;
    }

    // (first, second) is an essential pair that we cannot prove reduces
    // to zero.
    break; 
  }
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

  if (mTrackEssentialPair) {
    ASSERT(mDidReducingSPair.size() == newGen);
    mDidReducingSPair.push_back(false);
    if (!mTri.empty()) {
      ASSERT(mEssentialFirst < newGen ||
        mEssentialSecond == static_cast<size_t>(-1));
      if (mEssentialPoly.get() == 0 ||
        R->monomialIsDivisibleBy
          (mEssentialPoly->getLeadMonomial(), GB->getLeadMonomial(newGen))) {
        if (mEssentialFirst < newGen)
          --mEssentialSecond;
        nextEssentialPair();
      }
#ifdef DEBUG
      else {
        ASSERT(GB->basis().divisor(mEssentialPoly->getLeadMonomial()) ==
          static_cast<size_t>(-1));
      }
#endif
    }
  }

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

void SPairHandler::write(std::ostream &out) const
{
  out << "-- spairs --" << std::endl;
  // We write out all the pairs, one per line, in the form: [i j signature]
  // We do each element of the heap:
  for (size_t i = 0; i<heap.size(); i++)
    heap[i]->write(R, out);
}
void SPairHandler::dump() const
{
  write(std::cerr);
}

size_t SPairHandler::getMemoryUse() const
{
  return
    mTri.getMemoryUse() +
    heap.capacity() * sizeof(heap.front()) +
    getKnownSyzygyBitsMemoryUse();
}

size_t SPairHandler::getKnownSyzygyBitsMemoryUse() const {
  return mKnownSyzygyTri.getMemoryUse();
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
