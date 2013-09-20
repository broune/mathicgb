#include <utility>
#include <algorithm>
#include <iterator>

// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "mathicgb/stdinc.h"
#include "mathicgb/Range.hpp"

#include <gtest/gtest.h>
#include <iterator>
#include <vector>
#include <list>

using namespace mgb;

namespace {
  template<class Range1, class Range2>
  void checkRange(const Range1& r1, const Range2& r2) {
    const auto d1 = std::distance(std::begin(r1), std::end(r1));
    const auto d2 = std::distance(std::begin(r2), std::end(r2));
    ASSERT_EQ(d1, d2);
    ASSERT_TRUE(std::equal(std::begin(r1), std::end(r1), std::begin(r2)));
  }
}

TEST(Range, Simple) {
  int intArray[] = {1, 5, 3, 7};
  const std::vector<int> intVector(std::begin(intArray), std::end(intArray));
  const std::list<int> intList(std::begin(intArray), std::end(intArray));
  std::vector<int> v;

  auto checkClear = [&]() {
    checkRange(intVector, v);
    v.clear();
  };

  // Let's try the pattern without using range just to make sure everything
  // is OK.
  for (const auto& i : intArray)
    v.push_back(i);
  checkClear();

  // *** check using intArray.
  for (const auto& i : range(intArray))
    v.push_back(i);
  checkClear();

  for (const auto& i : range(std::begin(intArray), std::end(intArray)))
    v.push_back(i);
  checkClear();

  // *** check using intVector
  for (const auto& i : range(intVector))
    v.push_back(i);
  checkClear();

  for (const auto& i : range(std::begin(intVector), std::end(intVector)))
    v.push_back(i);
  checkClear();

  // *** check using intList
  for (const auto& i : range(intList))
    v.push_back(i);
  checkClear();

  for (const auto& i : range(std::begin(intList), std::end(intList)))
    v.push_back(i);
  checkClear();

  // *** check using iterated range
  for (const auto& i : range(range(std::begin(intArray), std::end(intArray))))
    v.push_back(i);
  checkClear();

  for (const auto& i : range(range(intVector)))
    v.push_back(i);
  checkClear();

  for (const auto& i : range(range(range(range(intList)))))
    v.push_back(i);
  checkClear();
}

TEST(Range, rangeToVector) {
  int intArray[] = {1, 5, 3, 7};
  const std::vector<int> intVector(std::begin(intArray), std::end(intArray));
  ASSERT_EQ(intVector, rangeToVector(std::begin(intArray), std::end(intArray)));
  ASSERT_EQ(intVector, rangeToVector(intArray));
  ASSERT_EQ(intVector, rangeToVector(range(intArray)));
}

TEST(Range, zip) {
  const std::string a[] = {"hello", "world"};
  const int b[] = {4, 2, 1, 0};
  std::ostringstream out;

  // Put a range() around b just to test something other than a built-in
  for (const auto& p : zip(a, range(b)))
    out << p.first << p.second << ' ';
  ASSERT_EQ("hello4 world2 ", out.str());

  // Now try the version on iterators and put the range around a instead.
  out.str("");
  for (
    const auto& p :
    zip(std::begin(range(a)), std::end(range(a)), std::begin(b), std::end(b))
  )
    out << p.first << p.second << ' ';
  ASSERT_EQ("hello4 world2 ", out.str());
}

TEST(Range, intRange) {
  const int int05[] = {0, 1, 2, 3, 4};
  ASSERT_EQ(rangeToVector(int05), rangeToVector(intRange(0, 5)));
  ASSERT_EQ(rangeToVector(int05), rangeToVector(intRange(5)));

  // gcc 4.7.3 won't parse "const signed char" here, but it will parse
  // it like this with the typedef.
  typedef signed char C;
  const C scharm2p5[] = {-2, -1, 0, 1, 2, 3, 4};
  ASSERT_EQ(
    rangeToVector(range(scharm2p5)),
    rangeToVector(intRange(C(-2), C(5)))
  );

  // Normally we should not dereference an end() iterator, but for the case
  // of intRange() it is OK.
  ASSERT_EQ(std::numeric_limits<size_t>::max(), *intRange().end());
  ASSERT_EQ(size_t(0), *intRange().begin());
}

TEST(Range, indexRange) {
  typedef std::pair<short, size_t> Pair;
  Pair indexed[] = {Pair(-2, 0), Pair(-1, 1), Pair(0, 2)};
  ASSERT_EQ(
    rangeToVector(indexed),
    rangeToVector(indexRange(intRange<short>(-2, 1)))
  );
}


TEST(Range, oppositePairRange) {
  std::string strs[] = {"hello", "world", "!"};
  const auto r = zip(strs, intRange(10));
  const auto opR = zip(intRange(3), strs);
  ASSERT_EQ(rangeToVector(r), rangeToVector(oppositePairRange(opR)));
  ASSERT_EQ(
    rangeToVector(opR),
    rangeToVector(oppositePairRange(std::begin(r), std::end(r)))
  );
}

TEST(Range, adjPairRange) {
  typedef std::pair<int, int> Pair;

  std::vector<Pair> none;
  ASSERT_EQ(rangeToVector(none), rangeToVector(adjPairRange(intRange(0))));
  ASSERT_EQ(rangeToVector(none), rangeToVector(adjPairRange(intRange(1))));

  Pair adj2[] = {Pair(0, 1)};
  ASSERT_EQ(rangeToVector(adj2), rangeToVector(adjPairRange(intRange(2))));

  Pair adj4[] = {Pair(0, 1), Pair(1, 2), Pair(2, 3)};
  ASSERT_EQ(rangeToVector(adj4), rangeToVector(adjPairRange(intRange(4))));
}
