// Copyright 2011 Michael E. Stillman

#ifndef _spair_handler_h_
#define _spair_handler_h_

// This class is designed for the Gao GB algorithm, or other signature based methods.
// The idea is to keep the size of the spair structures as small as possible

// Externally, an spair is (signature, integer).

#include "PairTriangle.hpp"
#include <vector>
#include "PolyRing.hpp"
#include "KoszulQueue.hpp"
#include <mathic.h>
#include <memtailor.h>

class Poly;
class MonomialTableArray;
class GroebnerBasis;
class FreeModuleOrder;
class Reducer;

class SPairHandler
{
public:
  SPairHandler(
    const PolyRing *R0,
    FreeModuleOrder *F0,
    const GroebnerBasis *GB0,
    MonomialTableArray *Hsyz0,
    Reducer* reducer,
    bool postponeKoszuls,
    bool useBaseDivisors,
    bool useSingularCriterionEarly,
    size_t queueType);
  ~SPairHandler();

  bool empty() const {return mTri.empty();}
  typedef std::vector<std::pair<size_t, size_t> > PairContainer;
  monomial popSignature(PairContainer& pairs);

  // fills in all the S-pairs with i.
  void newPairs(size_t i);

  // Inform the S-pair handler that there is a new syzygy signature in play.
  void newSyzygy(const_monomial sig);

  struct Stats {
    size_t comparisons; // comparisons not in construction (?)
    size_t precomparisons; // comparisons in spair construction (?)
    unsigned long long spairsConstructed; // all spairs
    unsigned long long spairsFinal; // spairs given to client
    unsigned long long nonregularSPairs; // spairs eliminated by being non-regular
    unsigned long long highBaseDivisorHits; // spairs eliminated by high base divisor
    unsigned long long lowBaseDivisorHits; // spairs eliminated by low base divisor
    unsigned long long hasHighBaseDivisor; // generators that have a high base divisor
    unsigned long long hasLowBaseDivisor; // generators that have a low base divisor
    unsigned long long syzygyModuleHits; // spairs eliminated by syzygy module
    unsigned long long earlyRelativelyPrimePairs;
    unsigned long long earlySingularCriterionPairs;
    unsigned long long queuedPairs; // number actually placed on spair triangle
    unsigned long long duplicateSignatures; // number of spairs removed due to duplicate signature

    Stats():
      comparisons(0),
      precomparisons(0),
      spairsConstructed(0),
      spairsFinal(0),
      nonregularSPairs(0),
      highBaseDivisorHits(0),
      lowBaseDivisorHits(0),
      hasHighBaseDivisor(0),
      hasLowBaseDivisor(0),
      syzygyModuleHits(0),
      earlyRelativelyPrimePairs(0),
      earlySingularCriterionPairs(0),
      queuedPairs(0),
      duplicateSignatures(0) {}
  };
  Stats getStats() const;

  size_t pairCount() const {return mTri.pairCount();}

  size_t getMemoryUse() const;
  size_t getKnownSyzygyBitsMemoryUse() const;

  // Informs the s-pair handler that the syzygy between gen1 and gen2
  // is a known syzygy.
  void setKnownSyzygy(size_t gen1, size_t gen2);
  void setKnownSyzygies(std::vector<std::pair<size_t, size_t> >& pairs);

  std::string name();

private:
  void makePreSPairs(size_t newGen);

  struct BaseDivisor { // a low ratio base divisor
    size_t baseDivisor; // the index of the generator that is the base divisor
    size_t ratioLessThan; // consider generators with ratio less than this
    monomial baseMonomial; // the monomial that has to divide to get a hit
  };
  typedef std::vector<BaseDivisor> BaseDivContainer;
  void setupBaseDivisors(
    BaseDivisor& divisor1,
    BaseDivisor& divisor2,
    size_t& highDivisorCmp,
    size_t newGenerator);
  
  const PolyRing *R;

  FreeModuleOrder *F;

  // if true, apply the early singular criterion
  bool const mUseSingularCriterionEarly;

  // true if low ratio base divisors are used to speed up S-pair elimination.
  const bool mUseBaseDivisors;

  // True if high ratio base divisors are used to speed up S-pair elimination.
  // The syzygy should have already been inserted into the syzygy module.
  const bool mUseHighBaseDivisors;

  // one entry for every s-pair, which is set to true if the
  // s-pair is known to be a syzygy. Only used if
  // mUseBaseDivisors is true.
  mathic::BitTriangle mKnownSyzygyTri;

  // From elsewhere
  MonomialTableArray *Hsyz; // we often modify this
  const GroebnerBasis *GB;
  Reducer* mReducer;
  const bool mPostponeKoszuls;

  class SigPairTriangle : public PairTriangle {
  public:
    SigPairTriangle(const GroebnerBasis& basis, size_t queueType);
  protected:
    virtual bool calculateOrderBy(size_t a, size_t b, monomial orderBy) const;
  private:
    const GroebnerBasis& mBasis;
  };
  SigPairTriangle mTri;

  mutable Stats mStats;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
