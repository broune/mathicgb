// Copyright 2011 Michael E. Stillman

#include "mathicgb/stdinc.h"

#include "mathicgb/Poly.hpp"
#include "mathicgb/Ideal.hpp"
#include "mathicgb/MTArray.hpp"
#include "mathicgb/MonTableNaive.hpp"
#include "mathicgb/io-util.hpp"
#include "mathicgb/PolyHeap.hpp"
#include "mathicgb/GroebnerBasis.hpp"
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

  std::unique_ptr<Ideal> I = idealParseFromString(idealA_fromStr_format);
  EXPECT_EQ("  -bc+ad\n  -b2+af\n  -bc2+a2e\n", toString(I.get()));
}

void testGB(int freeModuleOrder,
            std::string idealStr,
            std::string sigBasisStr,
            std::string syzygiesStr,
            std::string initialIdealStr,
            size_t nonSingularReductions)
{
  // Put the contents of pict.out into allPairsTest as a string. This
  // works because pict.out does not have any commas and we do not
  // care about whitespace. pict.out contains a set of tests such that
  // all pairs of parameters are covered by at least one test. See
  // pict.in for details.
#define MATHICGB_ESCAPE_MULTILINE_STRING(str) #str
  char const allPairsTests[] = MATHICGB_ESCAPE_MULTILINE_STRING(
spairQueue	reducerType	divLookup	monTable	buchberger	postponeKoszul	useBaseDivisors	autoTailReduce	autoTopReduce	preferSparseReducers	useSingularCriterionEarly
2	12	1	1	0	0	0	0	0	1	0
1	6	2	2	0	1	1	0	0	0	1
0	6	2	0	1	0	0	1	1	0	0
1	7	1	2	1	0	0	1	1	1	0
0	15	1	0	0	1	1	0	0	1	1
2	8	2	1	0	0	1	0	0	0	1
2	4	1	1	1	0	0	0	1	0	0
3	19	2	1	1	0	0	1	0	1	0
1	11	1	1	1	0	0	1	0	1	0
0	14	1	1	1	0	0	1	0	0	0
1	9	2	0	1	0	0	1	1	0	0
1	15	2	1	1	0	0	1	1	0	0
2	22	1	0	1	0	0	1	0	1	0
3	14	1	0	0	1	1	0	0	1	1
2	6	1	1	1	0	0	1	1	1	0
0	10	2	2	1	0	0	0	1	1	0
2	7	2	1	0	1	0	0	0	0	1
2	14	2	2	1	0	0	1	1	1	0
3	20	1	2	1	0	0	0	1	0	0
0	8	1	2	1	0	0	1	1	1	0
0	1	1	1	0	1	1	0	0	1	0
3	22	2	1	0	1	1	0	0	0	1
1	10	1	1	0	1	1	0	0	0	1
2	11	2	0	0	1	1	0	0	0	1
0	25	2	0	1	0	0	0	1	0	0
2	15	1	2	0	0	1	0	0	0	1
2	21	2	0	0	1	0	0	0	1	1
3	2	2	1	0	1	0	0	0	0	0
2	5	1	0	0	1	1	0	0	0	1
1	12	2	0	0	1	1	0	0	0	1
0	19	1	2	0	1	1	0	0	0	1
3	17	1	0	0	1	1	0	0	0	0
1	14	2	2	0	0	0	0	0	1	1
3	1	2	0	0	0	0	0	0	0	1
2	0	1	1	1	0	0	1	0	1	0
2	20	2	0	0	1	1	0	0	1	1
0	23	1	2	1	0	0	1	0	1	0
0	0	2	0	0	1	1	0	0	0	1
2	2	1	2	0	0	1	0	0	1	1
2	16	2	0	1	0	0	0	1	0	0
0	20	2	1	1	0	0	1	0	0	0
1	3	2	2	0	0	1	0	0	1	1
0	24	2	0	1	0	0	1	0	0	0
3	7	2	0	0	1	1	0	0	1	1
1	0	1	2	0	1	0	0	0	0	1
2	23	2	1	0	1	1	0	0	0	1
1	18	2	1	0	1	0	0	0	1	1
1	20	1	1	0	0	0	0	0	0	1
2	25	1	1	0	1	1	0	0	1	1
3	9	1	1	0	1	1	0	0	1	1
1	23	2	0	0	1	1	0	0	0	1
2	17	2	2	0	0	0	0	0	1	1
0	21	1	2	1	0	0	1	1	0	0
3	3	1	1	1	0	0	1	1	0	0
2	1	1	2	1	0	0	1	1	0	0
1	5	2	2	1	0	0	1	1	1	0
3	21	1	1	1	0	0	1	0	1	0
3	6	1	0	0	0	0	0	0	1	1
0	9	2	2	1	0	0	0	1	1	0
3	11	2	2	0	0	0	0	0	1	1
3	10	1	0	1	0	0	1	0	1	0
2	19	1	0	1	0	0	1	1	1	0
0	12	2	2	1	0	0	1	1	1	0
0	2	2	0	0	1	0	0	0	0	1
1	25	1	2	1	0	0	1	1	0	0
1	16	1	2	0	1	1	0	0	1	1
1	4	2	2	0	1	1	0	0	1	1
3	25	1	2	0	1	1	0	0	0	1
1	8	1	0	1	0	0	1	0	0	0
2	18	1	0	1	0	0	1	1	0	0
0	4	2	0	1	0	0	1	1	0	0
0	18	2	2	0	1	1	0	0	1	1
1	17	2	1	0	0	1	0	0	1	1
1	24	1	2	0	1	1	0	0	1	1
3	24	2	1	1	0	0	0	1	0	0
1	19	1	2	0	0	1	0	0	0	1
1	21	2	2	0	1	1	0	0	1	1
3	16	2	1	0	1	0	0	0	1	1
0	22	2	2	0	1	0	0	0	1	1
3	0	2	2	0	1	1	0	0	0	1
2	9	2	2	1	0	0	1	0	1	0
3	18	1	1	0	1	0	0	0	1	1
3	15	1	2	0	1	0	0	0	0	1
1	1	1	1	0	1	0	0	0	1	1
2	24	1	1	0	0	1	0	0	0	1
2	3	2	0	0	1	1	0	0	0	1
0	11	2	2	1	0	0	1	1	1	0
3	13	2	1	1	0	0	1	0	0	0
3	12	1	1	1	0	0	1	1	1	0
0	5	2	1	1	0	0	1	0	0	0
1	22	1	0	1	0	0	0	1	0	0
2	13	1	2	0	1	1	0	0	1	1
0	16	2	0	0	0	0	0	0	0	1
0	17	2	2	0	0	0	0	0	0	1
1	2	1	1	1	0	0	1	1	0	0
3	23	2	0	0	1	1	0	0	0	1
0	7	1	0	0	1	0	0	0	1	1
1	13	1	0	0	1	1	0	0	0	1
2	10	1	0	0	1	0	0	0	1	1
0	13	1	0	1	0	0	0	1	0	0
3	4	1	0	0	0	0	0	0	1	1
2	23	1	0	1	0	0	1	1	1	0
3	5	1	0	0	1	0	0	0	0	1
0	3	1	2	0	1	0	0	0	1	1
3	16	1	1	1	0	0	1	1	0	0
0	17	1	0	1	0	0	1	1	1	0
3	8	1	2	0	1	1	0	0	0	1
2	0	2	0	1	0	0	0	1	0	0
);
  std::istringstream tests(allPairsTests);
  // skip the initial line with the parameter names.
  {
      char const* params[] = {
        "spairQueue", "reducerType", "divLookup", "monTable",
        "buchberger", "postponeKoszul", "useBaseDivisors", "autoTailReduce",
        "autoTopReduce", "preferSparseReducers", "useSingularCriterionEarly"};

    std::string paramName;
    size_t const paramCount = sizeof(params) / sizeof(*params);
    for (size_t i = 0; i < paramCount; ++i) {
      tests >> paramName;
      // This assert will fire if you changed the order of the
      // parameters, renamed a parameter, removed a parameter or added
      // a parameter. Unless all you did was to rename a parameter,
      // don't just update the params array that the assert is based
      // on - you also need to update the code below that parses the
      // pict output because it depends on the order of the
      // parameters.
      MATHICGB_ASSERT(paramName == params[i]);
    }
  }

  while (true) {
    // parse a line of the pict file

    int spairQueue;
    tests >> spairQueue;
    if (!tests)
      break; // no more tests
    MATHICGB_ASSERT(0 <= spairQueue && spairQueue <= 3);

    int reducerType;
    tests >> reducerType;
    MATHICGB_ASSERT(0 <= reducerType && reducerType <= 30);

    int divLookup;
    tests >> divLookup;
    MATHICGB_ASSERT(1 <= divLookup && divLookup <= 2);

    int monTable;
    tests >> monTable;
    MATHICGB_ASSERT(0 <= monTable && monTable <= 2);
    
    int buchberger;
    tests >> buchberger;
    MATHICGB_ASSERT(0 <= buchberger && buchberger <= 1);

    int postponeKoszul;
    tests >> postponeKoszul;
    MATHICGB_ASSERT(0 <= postponeKoszul && postponeKoszul <= 1);

    int useBaseDivisors;
    tests >> useBaseDivisors;
    MATHICGB_ASSERT(0 <= useBaseDivisors && useBaseDivisors <= 1);

    int autoTailReduce;
    tests >> autoTailReduce;
    MATHICGB_ASSERT(0 <= autoTailReduce && autoTailReduce <= 1);

    int autoTopReduce;
    tests >> autoTopReduce;
    MATHICGB_ASSERT(0 <= autoTopReduce && autoTopReduce <= 1);

    int preferSparseReducers;
    tests >> preferSparseReducers;
    MATHICGB_ASSERT(0 <= preferSparseReducers && preferSparseReducers <= 1);

    int useSingularCriterionEarly;
    tests >> useSingularCriterionEarly;
    MATHICGB_ASSERT(0 <= useSingularCriterionEarly);
    MATHICGB_ASSERT(useSingularCriterionEarly <= 1);

    // Rule out combinations of parameter values that do not make sense.
    // These are asserts because pict should have already removed these
    // combinations.
    MATHICGB_ASSERT(buchberger || !autoTopReduce);
    MATHICGB_ASSERT(buchberger || !autoTailReduce);
    MATHICGB_ASSERT(!buchberger || !postponeKoszul);
    MATHICGB_ASSERT(!buchberger || !useBaseDivisors);
    MATHICGB_ASSERT(!buchberger || !useSingularCriterionEarly);

    // check that we have a valid reducer type
    Reducer::ReducerType red = Reducer::ReducerType(reducerType);
    MATHICGB_ASSERT(static_cast<int>(red) == reducerType);
    std::unique_ptr<Ideal> I(idealParseFromString(idealStr));
    MATHICGB_ASSERT
      (Reducer::makeReducerNullOnUnknown(red, I->ring()).get() != 0);

    if (buchberger) {
      BuchbergerAlg alg(
        *I, freeModuleOrder, Reducer::reducerType(reducerType), divLookup, preferSparseReducers, spairQueue);
      alg.setUseAutoTopReduction(autoTopReduce);
      alg.setUseAutoTailReduction(autoTailReduce);
      alg.computeGrobnerBasis();
      std::unique_ptr<Ideal> initialIdeal =
        alg.basis().initialIdeal();
      EXPECT_EQ(initialIdealStr, toString(initialIdeal.get()))
        << reducerType << ' ' << divLookup << ' '
        << monTable << ' ' << postponeKoszul << ' ' << useBaseDivisors;
    } else {
      SignatureGB basis
        (*I, freeModuleOrder, Reducer::reducerType(reducerType),
          divLookup, monTable, postponeKoszul, useBaseDivisors, preferSparseReducers, useSingularCriterionEarly, spairQueue);
      basis.computeGrobnerBasis();
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
