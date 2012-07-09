// Copyright 2011 Michael E. Stillman

#include "mathicgb/stdinc.h"

#include "mathicgb/Poly.hpp"
#include "mathicgb/Ideal.hpp"
#include "mathicgb/MTArray.hpp"
#include "mathicgb/MonTableNaive.hpp"
#include "mathicgb/io-util.hpp"
#include "mathicgb/PolyHeap.hpp"
#include "mathicgb/GroebnerBasis.hpp"
#include "mathicgb/SPairHandler.hpp"
#include "mathicgb/SignatureGB.hpp"
#include "mathicgb/BuchbergerAlg.hpp"
#include "test/ideals.hpp"

#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <memory>
#include <gtest/gtest.h>

extern int tracingLevel;

TEST(IO, ideal) {
  const char* idealA_fromStr_format = 
"32003 6 \
1 1 1 1 1 1 1 \
3 \
-bc+ad \
-b2+af \
-bc2+a2e \
";

  std::auto_ptr<Ideal> I = idealParseFromString(idealA_fromStr_format);
  EXPECT_EQ("  -bc+ad\n  -b2+af\n  -bc2+a2e\n", toString(I.get()));
}

void testGB(int freeModuleOrder,
            std::string idealStr,
            std::string sigBasisStr,
            std::string syzygiesStr,
            std::string initialIdealStr,
            size_t nonSingularReductions)
{
  for (int spairQueue = 0; spairQueue <= 3; ++spairQueue)
  for (int reducerType = 0; reducerType <= 30; ++reducerType)
  for (int divLookup = 1; divLookup <= 2; ++divLookup)
  for (int monTable = 0; monTable <= 2; ++monTable)
  for (int signatureBasis = 0; signatureBasis <= 1; ++signatureBasis)
  for (int buchberger = 0; buchberger <= 1; ++buchberger)
  for (int postponeKoszul = 0; postponeKoszul <= 1; ++postponeKoszul)
  for (int useBaseDivisors = 0; useBaseDivisors <= 1; ++useBaseDivisors)
  for (int autoTailReduce = 0; autoTailReduce <= 1; ++autoTailReduce)
  for (int autoTopReduce = 0; autoTopReduce <= 1; ++autoTopReduce)
  for (int preferSparseReducers = 0; preferSparseReducers <= 1;
    ++preferSparseReducers)
  for (int useSingularCriterionEarly = 0; useSingularCriterionEarly <= 1;
    ++useSingularCriterionEarly)
  {
    //std::cout << reducerType << ' ' << divLookup << ' ' << monTable << ' ' << signatureBasis << ' ' << buchberger << ' ' << postponeKoszul << ' ' << useBaseDivisors << ' ' << autoTailReduce << ' ' << autoTopReduce << ' ' << preferSparseReducers << std::endl;
    if (!buchberger && (autoTopReduce || autoTailReduce))
      continue;
    if (buchberger && (postponeKoszul || useBaseDivisors || signatureBasis || useSingularCriterionEarly))
      continue;

    Reducer::ReducerType red = Reducer::ReducerType(reducerType);
    if ((buchberger && signatureBasis) || static_cast<int>(red) != reducerType)
      continue;
    std::auto_ptr<Ideal> I(idealParseFromString(idealStr));
    if (Reducer::makeReducerNullOnUnknown(red, I->ring()).get() == 0)
      continue;

    if (buchberger) {
      ASSERT(!signatureBasis);
      BuchbergerAlg alg(
        *I, freeModuleOrder, Reducer::reducerType(reducerType), divLookup, preferSparseReducers, spairQueue);
      alg.setUseAutoTopReduction(autoTopReduce);
      alg.setUseAutoTailReduction(autoTailReduce);
      alg.computeGrobnerBasis();
      std::auto_ptr<Ideal> initialIdeal =
        alg.basis().initialIdeal();
      EXPECT_EQ(initialIdealStr, toString(initialIdeal.get()))
        << reducerType << ' ' << divLookup << ' '
        << monTable << ' ' << postponeKoszul << ' ' << useBaseDivisors;
    } else {
      SignatureGB basis
        (*I, freeModuleOrder, Reducer::reducerType(reducerType),
          divLookup, monTable, postponeKoszul, useBaseDivisors, preferSparseReducers, useSingularCriterionEarly, spairQueue);
      basis.setComputeSignatureBasis(signatureBasis);
      basis.computeGrobnerBasis();
      if (!signatureBasis) {
        std::auto_ptr<Ideal> initialIdeal =
          basis.getGB()->basis().initialIdeal();
        EXPECT_EQ(initialIdealStr, toString(initialIdeal.get()))
          << reducerType << ' ' << divLookup << ' '
          << monTable << ' ' << postponeKoszul << ' ' << useBaseDivisors;
      } else {
        EXPECT_EQ(sigBasisStr, toString(basis.getGB(), 1))
          << reducerType << ' ' << divLookup << ' '
          << monTable << ' ' << ' ' << postponeKoszul << ' '
          << useBaseDivisors;
        EXPECT_EQ(syzygiesStr, toString(basis.getSyzTable()))
          << reducerType << ' ' << divLookup << ' '
          << monTable << ' ' << ' ' << postponeKoszul << ' '
          << useBaseDivisors;
        EXPECT_EQ(nonSingularReductions, basis.getSigReductionCount() - basis.getSingularReductionCount())
          << reducerType << ' ' << divLookup << ' '
          << monTable << ' ' << ' ' << postponeKoszul << ' '
          << useBaseDivisors;
      }
    }
  }
}

extern int tracingLevel;
TEST(GB, small) {
  testGB(0, idealSmall, idealSmallBasis, idealSmallSyzygies, idealSmallInitial, 7);
}

TEST(GB, liu_0_1) {
  testGB(1, liu_ideal, liu_gb_strat0_free1,
    liu_syzygies_strat0_free1, liu_initial_strat0_free1, 13);
}

TEST(GB, weispfennig97_0_4) {
  testGB(4, weispfennig97_ideal, weispfennig97_gb_strat0_free4,
         weispfennig97_syzygies_strat0_free4, weispfennig97_initial_strat0_free4, 31);
}

TEST(GB, weispfennig97_0_5) {
  testGB(5, weispfennig97_ideal, weispfennig97_gb_strat0_free5,
         weispfennig97_syzygies_strat0_free5, weispfennig97_initial_strat0_free5, 27);
}

TEST(GB, gerdt93_0_1) {
  testGB(1, gerdt93_ideal, gerdt93_gb_strat0_free1,
         gerdt93_syzygies_strat0_free1, gerdt93_initial_strat0_free1, 9);
}

TEST(GB, gerdt93_0_2) {
  testGB(2, gerdt93_ideal, gerdt93_gb_strat0_free2,
         gerdt93_syzygies_strat0_free2, gerdt93_initial_strat0_free2, 7);
}

TEST(GB, gerdt93_0_3) {
  testGB(3, gerdt93_ideal, gerdt93_gb_strat0_free3,
         gerdt93_syzygies_strat0_free3, gerdt93_initial_strat0_free3, 9);
}

TEST(GB, gerdt93_0_4) {
  testGB(4, gerdt93_ideal, gerdt93_gb_strat0_free4,
         gerdt93_syzygies_strat0_free4, gerdt93_initial_strat0_free4, 7);
}
TEST(GB, gerdt93_0_5) {
  testGB(5, gerdt93_ideal, gerdt93_gb_strat0_free5,
         gerdt93_syzygies_strat0_free5, gerdt93_initial_strat0_free5, 7);
}


// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
