// Copyright 2011 Michael E. Stillman

#ifndef _sig_gb_h_
#define _sig_gb_h_

#include "PolyRing.hpp"
#include "MTArray.hpp"
#include "GroebnerBasis.hpp"
#include "FreeModuleOrder.hpp"
#include "SigSPairs.hpp"
#include "Reducer.hpp"
#include "KoszulQueue.hpp"
#include "SPairs.hpp"
#include "MonoProcessor.hpp"
#include <map>

class SigSPairs;
class DivisorLookup;

class SignatureGB {
public:
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::MonoVector MonoVector;

  SignatureGB(
    Basis&& basis,
    FreeModuleOrderType typ,
    Reducer::ReducerType reductiontyp,
    int divlookup_type,
    int montable_type,
    bool postponeKoszul,
    bool useBaseDivisors,
    bool preferSparseReducers,
    bool useSingularCriterionEarly,
    size_t queueType);

  void computeGrobnerBasis();

  // How many S-pairs were not eliminated before reduction.
  unsigned long long getSigReductionCount() const;

  // How many reductions were singular
  unsigned long long getSingularReductionCount() const;

  GroebnerBasis* getGB() { return GB.get(); }
  MonomialTableArray* getSyzTable() { return mProcessor->processingNeeded() ? Hsyz2.get() : Hsyz.get(); }
  SigSPairs* getSigSPairs() { return SP.get(); }

  size_t getMemoryUse() const;
  void displayStats(std::ostream& out) const;
  void displayPaperStats(std::ostream& out) const;
  void displayMemoryUse(std::ostream& out) const;
  void displaySomeStats(std::ostream& out) const;

  void setBreakAfter(unsigned int elements) {
    mBreakAfter = elements;
  }

  void setPrintInterval(unsigned int reductions) {
    mPrintInterval = reductions;
  }

  const Monoid& monoid() const {return R->monoid();}

private:
  unsigned int mBreakAfter;
  unsigned int mPrintInterval;




  bool processSPair(monomial sig, const SigSPairs::PairContainer& pairs);
  bool step();

  const PolyRing *R;
  std::unique_ptr<FreeModuleOrder> F;

  


  bool const mPostponeKoszul;

  // Currently we use either both criteria (high and loow) or neither.
  bool const mUseBaseDivisors;

  SigSPairs::PairContainer mSpairTmp; // use only for getting S-pairs

  // stats //////////
  size_t stats_sPairSignaturesDone; // distinct S-pair signatures done
  size_t stats_sPairsDone; // total S-pairs done
  size_t stats_koszulEliminated; // S-pairs eliminated due to Koszul queue
  size_t stats_SignatureCriterionLate; // # spairs removed due to being a syz signature

  // S-pairs eliminated due to relatively prime criterion
  size_t stats_relativelyPrimeEliminated;

  size_t stats_pairsReduced; // # spairs actually sent for reduction

  mic::Timer mTimer;
  double stats_nsecs;

  std::unique_ptr<GroebnerBasis> GB;
  std::unique_ptr<KoszulQueue> mKoszuls;
  std::unique_ptr<MonomialTableArray> Hsyz;
  std::unique_ptr<MonomialTableArray> Hsyz2;
  std::unique_ptr<Reducer> reducer;
  std::unique_ptr<SigSPairs> SP;
  std::unique_ptr<MonoProcessor<Monoid>> mProcessor;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
