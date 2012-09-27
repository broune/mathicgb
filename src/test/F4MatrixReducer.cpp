#include "mathicgb/stdinc.h"

#include "mathicgb/F4MatrixReducer.hpp"
#include "mathicgb/SparseMatrix.hpp"
#include "mathicgb/QuadMatrix.hpp"
#include "mathicgb/io-util.hpp"
#include "mathicgb/Poly.hpp"
#include <gtest/gtest.h>
#include <sstream>

TEST(F4MatrixReducer, Reduce) {
  auto ring = ringFromString("101 6 1\n1 1 1 1 1 1");

  QuadMatrix m;
  m.ring = ring.get();

  Poly p(ring.get());
  std::istringstream in("a1+a2+a3+a4+b1+b2+b3+b4+b5");
  p.parse(in);
  size_t count = 0;
  for (Poly::iterator it = p.begin(); it != p.end(); ++it) {
    monomial mono = it.getMonomial();
    if (count < 4)
      m.leftColumnMonomials.push_back(mono);
    else
      m.rightColumnMonomials.push_back(mono);
    ++count;
  }

  // top left
  m.topLeft.clear(4);
  m.topLeft.appendEntry(0, 1);
  m.topLeft.appendEntry(1, 2);
  m.topLeft.appendEntry(3, 3);
  m.topLeft.rowDone();
  m.topLeft.appendEntry(1, 1);
  m.topLeft.appendEntry(2, 3);
  m.topLeft.rowDone();
  m.topLeft.appendEntry(2, 1);
  m.topLeft.appendEntry(3, 7);
  m.topLeft.rowDone();
  m.topLeft.appendEntry(3, 1);
  m.topLeft.rowDone();

  // top right
  m.topRight.clear(5);
  m.topRight.appendEntry(2,8);
  m.topRight.rowDone();
  m.topRight.appendEntry(3,9);
  m.topRight.rowDone();
  m.topRight.appendEntry(4,10);
  m.topRight.rowDone();
  m.topRight.rowDone();

  // bottom left
  m.bottomLeft.clear(4);
  m.bottomLeft.appendEntry(0, 1);
  m.bottomLeft.appendEntry(1, 1);
  m.bottomLeft.appendEntry(3, 24);
  m.bottomLeft.rowDone();
  m.bottomLeft.appendEntry(1, 9);
  m.bottomLeft.rowDone();

  // bottom right
  m.bottomRight.clear(5);
  m.bottomRight.appendEntry(0, 1);
  m.bottomRight.appendEntry(2, 12);
  m.bottomRight.appendEntry(3, 13);
  m.bottomRight.appendEntry(4, 41);
  m.bottomRight.rowDone();
  m.bottomRight.appendEntry(1, 2);
  m.bottomRight.appendEntry(3, 11);
  m.bottomRight.rowDone();

  MATHICGB_ASSERT(m.debugAssertValid());
  const char* origStr = 
    "Left columns: a a2 a3 a4\n"
    "Right columns: b b2 b3 b4 b5\n"
    "0: 0#1 1#2 3#3  | 0: 2#8               \n"
    "1: 1#1 2#3      | 1: 3#9               \n"
    "2: 2#1 3#7      | 2: 4#10              \n"
    "3: 3#1          | 3:                   \n"
    "                |                      \n"
    "0: 0#1 1#1 3#24 | 0: 0#1 2#12 3#13 4#41\n"
    "1: 1#9          | 1: 1#2 3#11          \n";
  ASSERT_EQ(origStr, m.toString());

  SparseMatrix reduced;
  F4MatrixReducer red;
  red.reduce(*ring, m, reduced);

  const char* redStr =
    "0: 0#1 2#4 3#22 4#11\n"
    "1: 1#1 3#66 4#34\n";
  ASSERT_EQ(redStr, reduced.toString());
}
