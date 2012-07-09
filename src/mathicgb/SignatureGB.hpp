// Copyright 2011 Michael E. Stillman

#ifndef _sig_gb_h_
#define _sig_gb_h_

#ifdef Win32
#include "StdAfx.h"
#endif

#include "PolyRing.hpp"
#include "MTArray.hpp"
#include "GroebnerBasis.hpp"
#include "FreeModuleOrder.hpp"
#include "SPairHandler.hpp"
#include "Reducer.hpp"
#include "KoszulQueue.hpp"
#include "SPairs.hpp"
#include <map>

class SPairHandler;
class DivisorLookup;

class SignatureGB {
public:
  SignatureGB(
    const Ideal& ideal,
    FreeModuleOrderType typ,
    Reducer::ReducerType reductiontyp,
    int divlookup_type,
    int montable_type,
    bool postponeKoszul,
    bool useBaseDivisors,
    bool preferSparseReducers,
    bool useSingularCriterionEarly,
    size_t queueType);
  ~SignatureGB();

  void computeGrobnerBasis();

  // How many S-pairs were not eliminated before reduction.
  unsigned long long getSigReductionCount() const;

  // How many reductions were singular
  unsigned long long getSingularReductionCount() const;

  GroebnerBasis* getGB() { return GB; }
  MonomialTableArray* getSyzTable() { return Hsyz; }
  SPairHandler* getSPairHandler() { return SP; }

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

  void setComputeSignatureBasis(bool value) {
    mComputeSignatureBasis = value;
  }

private:
  unsigned int mBreakAfter;
  unsigned int mPrintInterval;
  bool mComputeSignatureBasis;

  bool processSPair(monomial sig, const SPairHandler::PairContainer& pairs);
  void step();

  const PolyRing *R;
  FreeModuleOrder *F;

  SPairHandler *SP;
  MonomialTableArray *Hsyz;
  GroebnerBasis *GB;
  KoszulQueue *mKoszuls;

  Reducer* reducer;

  bool const mPostponeKoszul;

  // Currently we use either both criteria (high and loow) or neither.
  bool const mUseBaseDivisors;

  SPairHandler::PairContainer mSpairTmp; // use only for getting S-pairs

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
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
