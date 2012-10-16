#include "mathicgb/stdinc.h"

#include "mathicgb/SparseMatrix.hpp"
#include "mathicgb/Poly.hpp"
#include "mathicgb/PolyRing.hpp"
#include "mathicgb/io-util.hpp"
#include <gtest/gtest.h>
#include <memory>

namespace {
  std::unique_ptr<Poly> parsePoly(const PolyRing& ring, std::string str) {
    auto p = make_unique<Poly>(ring);
    std::istringstream in(str);
    p->parse(in);
    return p;
  }
}

TEST(SparseMatrix, NoRows) {
  SparseMatrix mat; // test a matrix with no rows
  ASSERT_EQ(0, mat.entryCount());
  ASSERT_EQ(0, mat.rowCount());
  ASSERT_EQ(0, mat.colCount());
  ASSERT_EQ("matrix with no rows\n", mat.toString()); 
}

TEST(SparseMatrix, Simple) {
  SparseMatrix mat(2000);
  mat.appendColumn();

  mat.appendEntry(5, 101);
  mat.rowDone();
  ASSERT_EQ(1, mat.entryCount());
  ASSERT_EQ(1, mat.rowCount());
  ASSERT_EQ(2001, mat.colCount());
  ASSERT_EQ(5, mat.leadCol(0));
  ASSERT_EQ(1, mat.entryCountInRow(0));
  ASSERT_EQ("0: 5#101\n", mat.toString()); 
  ASSERT_FALSE(mat.emptyRow(0));

  mat.rowDone(); // add a row with no entries
  ASSERT_EQ(1, mat.entryCount());
  ASSERT_EQ(2, mat.rowCount());
  ASSERT_EQ(2001, mat.colCount());
  ASSERT_EQ(5, mat.leadCol(0));
  ASSERT_EQ(0, mat.entryCountInRow(1));
  ASSERT_EQ("0: 5#101\n1:\n", mat.toString()); 
  ASSERT_TRUE(mat.emptyRow(1));

  mat.appendEntry(5, 102);
  mat.appendColumn();
  mat.appendEntry(2001, 0); // scalar zero
  mat.rowDone(); // add a row with two entries
  ASSERT_EQ(3, mat.entryCount());
  ASSERT_EQ(3, mat.rowCount());
  ASSERT_EQ(2002, mat.colCount());
  ASSERT_EQ(5, mat.leadCol(2));
  ASSERT_EQ(2, mat.entryCountInRow(2));
  ASSERT_EQ("0: 5#101\n1:\n2: 5#102 2001#0\n", mat.toString()); 
  ASSERT_FALSE(mat.emptyRow(2));
}

TEST(SparseMatrix, toRow) {
  auto ring = ringFromString("32003 6 1\n1 1 1 1 1 1");
  auto polyForMonomials = parsePoly(*ring, "a5+a4+a3+a2+a1+a0");
  std::vector<monomial> monomials;
  for (auto it = polyForMonomials->begin(); it != polyForMonomials->end(); ++it)
    monomials.push_back(it.getMonomial());

  SparseMatrix mat(5);
  mat.clear(6);
  mat.rowDone();
  mat.appendEntry(0,10);
  mat.rowDone();
  mat.appendEntry(2,20);
  mat.appendEntry(3,0);
  mat.appendEntry(4,40);
  mat.rowDone();

  Poly p(*ring);
  mat.rowToPolynomial(0, monomials, p);
  ASSERT_EQ(*parsePoly(*ring, "0"), p);
  mat.rowToPolynomial(1, monomials, p);
  ASSERT_EQ(*parsePoly(*ring, "10a5"), p);
  mat.rowToPolynomial(2, monomials, p);
  ASSERT_EQ(*parsePoly(*ring, "20a3+40a1"), p);
}
