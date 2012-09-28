#include "stdinc.h"
#include "BuchbergerAlg.hpp"
#include "Ideal.hpp"

#include <iostream>

extern int tracingLevel;

BuchbergerAlg::BuchbergerAlg(
  const Ideal& ideal,
  FreeModuleOrderType orderType,
  Reducer::ReducerType reducerType,
  int divisorLookupType,
  bool preferSparseReducers,
  size_t queueType
):
  mBreakAfter(0),
  mPrintInterval(0),
  mSPairGroupSize(0),
  mUseAutoTopReduction(true),
  mUseAutoTailReduction(false),
  mRing(*ideal.getPolyRing()),
  mOrder(FreeModuleOrder::makeOrder(orderType, &ideal)),
  mReducer(Reducer::makeReducer(reducerType, *ideal.getPolyRing())),
  mBasis(mRing, *mOrder, DivisorLookup::makeFactory(
    *ideal.getPolyRing(),
    divisorLookupType)->create(preferSparseReducers, true)
  ),
  mSPairs(mBasis, queueType),
  mSPolyReductionCount(0)
{
  // Reduce and insert the generators of the ideal into the starting basis
  size_t const idealSize = ideal.size();
  for (size_t gen = 0; gen != idealSize; ++gen)
    insertReducedPoly(mReducer->classicReduce(*ideal.getPoly(gen), mBasis));
}

void BuchbergerAlg::insertReducedPoly(
  std::unique_ptr<Poly> polyToInsert
) {
  ASSERT(polyToInsert.get() != 0);
  if (polyToInsert->isZero())
    return;
  ASSERT(mBasis.divisor(polyToInsert->getLeadMonomial()) ==
    static_cast<size_t>(-1));

  if (tracingLevel > 20) {
    std::cerr << "inserting basis element " << mBasis.size() << ": ";
    if (tracingLevel > 100) {
      polyToInsert->display(std::cerr);
      std::cerr << std::endl;
    } else {
      mRing.printMonomialFrobbyM2Format
        (std::cerr, polyToInsert->getLeadMonomial());
      if (polyToInsert->nTerms() > 1)
        std::cerr << " + [...]";
      std::cerr << std::endl;
    }
  }

  if (!mUseAutoTopReduction) {
    size_t const newGen = mBasis.size();
    mBasis.insert(std::move(polyToInsert));
    mSPairs.addPairs(newGen);
    return;
  }

  std::vector<size_t> toRetireAndReduce;
  std::vector<Poly*> toReduce;

  try {
    do {
      // reduce polynomial and insert into basis
      {
        std::unique_ptr<Poly> reduced;
        if (toReduce.empty()) // if first iteration
          reduced = std::move(polyToInsert);
        else {
          reduced = mReducer->classicReduce(*toReduce.back(), mBasis);
          if (tracingLevel > 20) {
            if (reduced->isZero()) {
              std::cerr << "auto-top-reduce cascade: "
                "basis element reduced to zero."
                << std::endl;
            } else {
              std::cerr << "auto-top-reduce cascade: "
                "inserting reduced poly with lead term "
                << std::endl;
              mRing.printMonomialFrobbyM2Format
              (std::cerr, reduced->getLeadMonomial());
              std::cerr << '\n';
            }
          }
          delete toReduce.back();
          toReduce.pop_back();
        }
        if (reduced->isZero())
          continue;
        reduced->makeMonic(); 
        mBasis.insert(std::move(reduced));
      }

      // form S-pairs and retire basis elements that become top reducible.
      const size_t newGen = mBasis.size() - 1;
      mSPairs.addPairsAssumeAutoReduce(newGen, toRetireAndReduce);
      for (std::vector<size_t>::const_iterator it = toRetireAndReduce.begin();
        it != toRetireAndReduce.end(); ++it) {
        toReduce.push_back(0); // allocate space in vector before .release()
        toReduce.back() = mBasis.retire(*it).release();
      }
      toRetireAndReduce.clear();
    } while (!toReduce.empty());
  } catch (...) {
    for (std::vector<Poly*>::iterator it = toReduce.begin();
      it != toReduce.end(); ++it)
      delete *it;
    throw;
  }
  ASSERT(toReduce.empty());
  ASSERT(toRetireAndReduce.empty());
}

void BuchbergerAlg::computeGrobnerBasis() {
  size_t counter = 0;
  mTimer.reset();
  mRing.resetCoefficientStats();
  
  if (mUseAutoTailReduction)
    autoTailReduce();

  while (!mSPairs.empty()) {
    step();
    if (mBreakAfter != 0 && mBasis.size() > mBreakAfter) {
      std::cerr
        << "Stopping Grobner basis computation due to reaching limit of "
        << mBreakAfter << " basis elements." << std::endl;
      break;
    }
    if (mPrintInterval != 0 && (++counter % mPrintInterval) == 0)
      printStats(std::cerr);
  }
  //printStats(std::cerr);
  //mReducer->dump();

  for (size_t i = 0; i < mBasis.size(); ++i)
    if (!mBasis.retired(i))
      mBasis.replaceSameLeadTerm
        (i, mReducer->classicTailReduce(mBasis.poly(i), mBasis));
}

void BuchbergerAlg::step() {
  ASSERT(!mSPairs.empty());
  if (tracingLevel > 30)
    std::cerr << "Determining next S-pair" << std::endl;

  if (mSPairGroupSize == 0) {
    std::pair<size_t, size_t> p = mSPairs.pop();
    if (p.first == static_cast<size_t>(-1)) {
      ASSERT(p.second == static_cast<size_t>(-1));
      return; // no more S-pairs
    }
    ASSERT(p.first != static_cast<size_t>(-1));
    ASSERT(p.second != static_cast<size_t>(-1));
    ASSERT(!mBasis.retired(p.first));
    ASSERT(!mBasis.retired(p.second));

    if (tracingLevel > 20) {
      std::cerr << "Reducing S-pair ("
                << p.first << ", "
                << p.second << ")" << std::endl;
    }
    std::unique_ptr<Poly> reduced
      (mReducer->classicReduceSPoly
       (mBasis.poly(p.first), mBasis.poly(p.second), mBasis));
    if (!reduced->isZero()) {
      insertReducedPoly(std::move(reduced));
      if (mUseAutoTailReduction)
        autoTailReduce();
    }
  } else {
    std::vector<std::pair<size_t, size_t> > spairGroup;
    for (unsigned int i = 0; i < mSPairGroupSize; ++i) {
      std::pair<size_t, size_t> p = mSPairs.pop();
      if (p.first == static_cast<size_t>(-1)) {
        ASSERT(p.second == static_cast<size_t>(-1));
        break; // no more S-pairs
      }
      ASSERT(p.first != static_cast<size_t>(-1));
      ASSERT(p.second != static_cast<size_t>(-1));
      ASSERT(!mBasis.retired(p.first));
      ASSERT(!mBasis.retired(p.second));

      spairGroup.push_back(p);
    }
    if (spairGroup.empty())
      return; // no more s-pairs
    std::vector<std::unique_ptr<Poly> > reduced;
    mReducer->classicReduceSPolyGroup(spairGroup, mBasis, reduced);
    
    for (auto it = reduced.begin(); it != reduced.end(); ++it) {
      auto p = std::move(*it);
      MATHICGB_ASSERT(!p->isZero());
      if (it != reduced.begin())
        p = mReducer->classicReduce(*p, mBasis);
      if (!p->isZero())
        insertReducedPoly(std::move(p));
    }
    if (mUseAutoTailReduction)
      autoTailReduce();
  }
}

void BuchbergerAlg::autoTailReduce() {
  ASSERT(mUseAutoTailReduction);

  for (size_t i = 0; i < mBasis.size(); ++i) {
    if (mBasis.retired(i))
      continue;
    if (mBasis.usedAsReducerCount(i) < 1000)
      continue;
    mBasis.replaceSameLeadTerm
      (i, mReducer->classicTailReduce(mBasis.poly(i), mBasis));
  }
}

size_t BuchbergerAlg::getMemoryUse() const {
  return
    mBasis.getMemoryUse() +
    mRing.getMonomialPool().getMemoryUse() +
    mReducer->getMemoryUse() +
    mSPairs.getMemoryUse();
}

void BuchbergerAlg::printStats(std::ostream& out) const {
  out << " term order:         " << mBasis.order().description() << '\n';
  out << " reduction type:     " << mReducer->description() << '\n';
  out << " divisor tab type:   " << mBasis.divisorLookup().getName() << '\n';
  out << " S-pair queue type:  " << mSPairs.name() << '\n';
  out << " total compute time: " << mTimer.getMilliseconds()/1000.0 << " seconds " << '\n';
  out << " S-pair group size:  " << mSPairGroupSize << '\n';

  const PolyRing::coefficientStats& cstats = mRing.getCoefficientStats();
  out << "n-coeff-add:         " << cstats.n_add << '\n';
  out << "n-coeff-addmult:     " << cstats.n_addmult << '\n';
  out << "n-coeff-mult:        " << cstats.n_mult << '\n';
  out << "n-coeff-recip:       " << cstats.n_recip << '\n';
  out << "n-coeff-divide:      " << cstats.n_divide << '\n';

  mic::ColumnPrinter pr;
  pr.addColumn(true, " ");
  pr.addColumn(false, " ");
  pr.addColumn(true, " ");

  std::ostream& name = pr[0];
  std::ostream& value = pr[1];
  std::ostream& extra = pr[2];

  const size_t basisSize = mBasis.size();
  const double mseconds = mTimer.getMilliseconds();
  const size_t pending = mSPairs.pairCount();

  name << "Time spent:\n";
  value << mTimer << '\n';
  extra << mic::ColumnPrinter::oneDecimal(mseconds / basisSize)
        << " ms per basis element\n";

  const double pendingRatio = static_cast<double>(pending) / basisSize;
  name << "Basis elements:\n";
  value << mic::ColumnPrinter::commafy(basisSize) << '\n';
  extra << mic::ColumnPrinter::oneDecimal(pendingRatio)
        << " Sp pend per basis ele\n";

  const size_t basisTermCount = mBasis.monomialCount();
  name << "Terms for basis:\n";
  value << mic::ColumnPrinter::commafy(basisTermCount) << '\n';
  extra << mic::ColumnPrinter::ratio(basisTermCount, basisSize)
        << " terms per basis ele\n";

  const size_t minLeadCount = mBasis.minimalLeadCount();
  name << "Minimum lead terms:\n";
  value << mic::ColumnPrinter::commafy(minLeadCount) << '\n';
  extra << mic::ColumnPrinter::percent(minLeadCount, basisSize)
        << " basis ele have min lead\n";

  const size_t lastMinLead = mBasis.maxIndexMinimalLead() + 1;
  const size_t timeSinceLastMinLead = basisSize - lastMinLead;
  name << "Index of last min lead:\n";
  value << mic::ColumnPrinter::commafy(lastMinLead) << '\n';
  extra << mic::ColumnPrinter::percent(timeSinceLastMinLead, basisSize)
        << " of basis added since then\n";

  const unsigned long long considered =
    mBasis.size() * (mBasis.size() - 1) / 2;
  name << "S-pairs considered:\n";
  value << mic::ColumnPrinter::commafy(considered) << '\n';
  extra << '\n';

  name << "S-pairs pending:\n";
  value << mic::ColumnPrinter::commafy(pending) << '\n';
  extra << mic::ColumnPrinter::percent(pending, considered)
        << " of considered\n";

  unsigned long long const reductions = sPolyReductionCount();
  name << "S-pairs reduced:\n";
  value << mic::ColumnPrinter::commafy(reductions) << '\n';
  extra << '\n';

  Reducer::Stats reducerStats = mReducer->sigStats();
  SPairs::Stats sPairStats = mSPairs.stats();

  unsigned long long const primeElim = sPairStats.relativelyPrimeHits;
  name << "Rel.prime sp eliminated:\n";
  value << mic::ColumnPrinter::commafy(primeElim) << '\n';
  extra << mic::ColumnPrinter::percent(primeElim, reductions)
        << " of late eliminations\n";

  const unsigned long long singularReductions =
    reducerStats.singularReductions;
  name << "Singular reductions:\n";
  value << mic::ColumnPrinter::commafy(singularReductions) << '\n';
  extra << mic::ColumnPrinter::percent(singularReductions, reductions)
        << " of reductions\n";

  const unsigned long long zeroReductions = reducerStats.zeroReductions;
  name << "Reductions to zero:\n";
  value << mic::ColumnPrinter::commafy(zeroReductions) << '\n';
  extra << mic::ColumnPrinter::percent(zeroReductions, reductions)
        << " of reductions\n";

  const unsigned long long newReductions =
    reductions - singularReductions - zeroReductions;
  name << "Reductions to new ele:\n";
  value << mic::ColumnPrinter::commafy(newReductions) << '\n';
  extra << mic::ColumnPrinter::percent(newReductions, reductions)
        << " of reductions\n";

  const unsigned long long redSteps = reducerStats.steps;
  const double stepsRatio =
    static_cast<double>(redSteps) / (reductions - singularReductions);
  name << "Sig reduction steps:\n";
  value << mic::ColumnPrinter::commafy(redSteps) << '\n';
  extra << mic::ColumnPrinter::oneDecimal(stepsRatio)
        << " steps per non-sing reduction\n";

  const unsigned long long longestReduction = reducerStats.maxSteps;
  name << "Longest sig reduction:\n";
  value << mic::ColumnPrinter::commafy(longestReduction) << '\n';
  extra << '\n';

  Reducer::Stats classicRedStats = mReducer->classicStats();
  const unsigned long long clReductions = classicRedStats.reductions;
  name << "Classic reductions:\n";
  value << mic::ColumnPrinter::commafy(clReductions) << '\n';
  extra << '\n';

  const unsigned long long clRedSteps =  classicRedStats.steps;
  const double clStepsRatio = static_cast<double>(clRedSteps) / clReductions;
  name << "Classic reduction steps:\n";
  value << mic::ColumnPrinter::commafy(clRedSteps) << '\n';
  extra << mic::ColumnPrinter::oneDecimal(clStepsRatio)
        << " steps per reduction\n";

  const unsigned long long clLongestReduction = classicRedStats.maxSteps;
  name << "Longest classic red:\n";
  value << mic::ColumnPrinter::commafy(clLongestReduction) << '\n';
  extra << '\n';

  //SPairs::Stats sPairStats = mSPairs.stats();  
  unsigned long long marginal = sPairStats.sPairsConsidered;

  unsigned long long const primeHits = sPairStats.relativelyPrimeHits;
  name << "Buchb relatively prime:\n";
  value << mic::ColumnPrinter::commafy(primeHits) << '\n';
  extra << mic::ColumnPrinter::percent(primeHits, marginal) <<
    " of S-pairs\n";
  marginal -= primeHits;

  unsigned long long const simpleHits = sPairStats.buchbergerLcmSimpleHits;
  name << "Buchb lcm simple hits:\n";
  value << mic::ColumnPrinter::commafy(simpleHits) << '\n';
  extra << mic::ColumnPrinter::percent(simpleHits, marginal) <<
    " of remaining S-pairs\n";
  marginal -= simpleHits;

  unsigned long long const simpleHitsLate = sPairStats.buchbergerLcmSimpleHitsLate;
  name << "Buchb late lcm simple hits:\n";
  value << mic::ColumnPrinter::commafy(simpleHitsLate) << '\n';
  extra << mic::ColumnPrinter::percent(simpleHitsLate, marginal) <<
    " of remaining S-pairs\n";
  marginal -= simpleHitsLate;

  unsigned long long const buchAdvHits =
    sPairStats.buchbergerLcmAdvancedHits;
  name << "Buchb lcm adv hits:\n";
  value << mic::ColumnPrinter::commafy(buchAdvHits) << '\n';
  extra << mic::ColumnPrinter::percent(buchAdvHits, marginal) <<
    " of remaining S-pairs\n";

  const unsigned long long buchCache = sPairStats.buchbergerLcmCacheHits;
  name << "Buchb lcm cache hits:\n";
  value << mic::ColumnPrinter::commafy(buchCache) << '\n';
  extra << mic::ColumnPrinter::percent(buchCache, simpleHits) <<
    " of simple hits\n";

  const unsigned long long buchCacheLate = sPairStats.buchbergerLcmCacheHitsLate;
  name << "Buchb late lcm cache hits:\n";
  value << mic::ColumnPrinter::commafy(buchCacheLate) << '\n';
  extra << mic::ColumnPrinter::percent(buchCacheLate, simpleHits) <<
    " of simple hits\n";

  out << "***** Classic Buchberger algorithm statistics *****\n"
      << pr << std::flush;
}

void BuchbergerAlg::printMemoryUse(std::ostream& out) const
{
  // Set up printer
  mic::ColumnPrinter pr;
  pr.addColumn();
  pr.addColumn(false);
  pr.addColumn(false);

  std::ostream& name = pr[0];
  std::ostream& value = pr[1];
  std::ostream& extra = pr[2];

  const size_t total = getMemoryUse();
  { // Grobner basis
    const size_t basisMem = mBasis.getMemoryUse();
    name << "Grobner basis:\n";
    value << mic::ColumnPrinter::bytesInUnit(basisMem) << '\n';
    extra << mic::ColumnPrinter::percent(basisMem, total) << '\n';
  }
  { // Spairs
    const size_t sPairMem = mSPairs.getMemoryUse();
    name << "S-pairs:\n";
    value << mic::ColumnPrinter::bytesInUnit(sPairMem) << '\n';
    extra << mic::ColumnPrinter::percent(sPairMem, total) << '\n';
  }
  { // Reducer
    const size_t reducerMem = mReducer->getMemoryUse();
    name << "Reducer:\n";
    value << mic::ColumnPrinter::bytesInUnit(reducerMem) << '\n';
    extra << mic::ColumnPrinter::percent(reducerMem, total) << '\n';
  }
  { // Signatures
    const size_t sigMem = mRing.getMonomialPool().getMemoryUse();
    name << "Monomials:\n";
    value << mic::ColumnPrinter::bytesInUnit(sigMem) << '\n';
    extra << mic::ColumnPrinter::percent(sigMem, total) << '\n';
  }
  // total
  name << "-------------\n";
  value << '\n';
  extra << '\n';

  name << "Memory used in total:\n";
  value << mic::ColumnPrinter::bytesInUnit(total) << "\n";
  extra << "\n";

  out << "*** Memory use by component ***\n" << pr << std::flush;
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
