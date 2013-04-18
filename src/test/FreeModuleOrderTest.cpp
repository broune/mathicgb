#include "mathicgb/stdinc.h"

#include "mathicgb/PolyRing.hpp"
#include "mathicgb/Ideal.hpp"
#include "mathicgb/FreeModuleOrder.hpp"
#include "mathicgb/io-util.hpp"
#include <gtest/gtest.h>
#include <algorithm>

void runTest(
  const char* idealStr,
  const char* signatureStr,
  const char* correctStr,
  int orderType
) {
  std::string line;

  std::unique_ptr<Ideal> ideal = idealParseFromString(idealStr);
  const PolyRing* ring = ideal->getPolyRing();

  std::vector<monomial> sigs;
  std::vector<PreSPair> pairs;
  {
    std::istringstream in(signatureStr);
    while (std::getline(in, line)) {
      sigs.push_back(monomialParseFromString(ring, line));
      PreSPair pair;
      pair.i = static_cast<BigIndex>(pairs.size());
      pair.signature = monomialParseFromString(ring, line);
      pairs.push_back(pair);
    }
  }

  std::vector<size_t> answer(pairs.size());
  {
    std::istringstream in(correctStr);
    for (size_t i = 0; i < answer.size(); ++i) {
      in >> answer[i];
      MATHICGB_ASSERT(!!in);
    }
  }

  MATHICGB_ASSERT(sigs.size() == pairs.size());
  MATHICGB_ASSERT(sigs.size() == pairs.size());

  std::unique_ptr<FreeModuleOrder> order
    (FreeModuleOrder::makeOrder(orderType, ideal.get()));
  order->sortAndScrambleSignatures(pairs);
  for (size_t i = 0; i < pairs.size(); ++i) {
    ring->freeMonomial(pairs[i].signature);
    pairs[i].signature = sigs[pairs[i].i];
  }
  sigs.clear();

  for (size_t i = 0; i < pairs.size(); ++i)
    ASSERT_EQ(pairs[i].i, answer[i]) << i << ' ' << orderType;
  for (size_t i = 0; i < pairs.size(); ++i) {
    const_monomial sigi = pairs[i].signature;
    ASSERT_EQ(order->signatureCompare(sigi, sigi), EQ)
      << i << ' ' << orderType;
    for (size_t j = 0; j < i; ++j) {
      const_monomial sigj = pairs[j].signature;
      ASSERT_EQ(order->signatureCompare(sigj, sigi), LT)
         << i << ' ' << j << ' ' << orderType;
      ASSERT_EQ(order->signatureCompare(sigi, sigj), GT)
         << i << ' ' << j << ' ' << orderType;
    }
  }
}

TEST(FreeModuleOrder, One) {
  const char* ideal = 
"32003 3 "
"1 1 1 1 "
"3 "
"a3b2c2 "
"a "
"a3c2 ";
  const char* sigs =
    "<1>\n" // 0
    "<0>\n" // 1
    "b10<0>\n" // 2
    "ac<0>\n" // 3
    "bc<0>\n" // 4
    "ac<1>\n" // 5
    "bc<1>\n" // 6
    "ab2c<2>\n"; // 7

  runTest(ideal, sigs, "0 1 6 4 5 3 7 2", 1);
  //runTest(ideal, sigs, "0 6 5 1 4 3 7 2", 2); 
  runTest(ideal, sigs, "0 6 5 1 7 4 3 2", 3); 
  //runTest(ideal, sigs, "0 6 5 1 4 3 7 2", 4); 
  //runTest(ideal, sigs, "0 6 5 1 4 7 3 2", 5); 
  runTest(ideal, sigs, "1 4 3 2 0 6 5 7", 6);
  runTest(ideal, sigs, "7 0 6 5 1 4 3 2", 7);
}
