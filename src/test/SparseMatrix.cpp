#include "mathicgb/stdinc.h"

#include "mathicgb/SparseMatrix.hpp"
#include <gtest/gtest.h>

TEST(SparseMatrix, NoRows) {
  SparseMatrix mat; // test a matrix with no rows
  ASSERT_EQ(0, mat.entryCount());
  ASSERT_EQ(0, mat.rowCount());
  ASSERT_EQ(0, mat.colCount());
  ASSERT_EQ("matrix with no rows\n", mat.toString()); 
}

TEST(SparseMatrix, Simple) {
  SparseMatrix mat;

  mat.appendEntry(5, 101);
  mat.rowDone();
  ASSERT_EQ(1, mat.entryCount());
  ASSERT_EQ(1, mat.rowCount());
  ASSERT_EQ(6, mat.colCount());
  ASSERT_EQ(5, mat.leadCol(0));
  ASSERT_EQ(1, mat.entryCountInRow(0));
  ASSERT_EQ("0: 5#101\n", mat.toString()); 
  ASSERT_FALSE(mat.emptyRow(0));

  mat.rowDone(); // add a row with no entries
  ASSERT_EQ(1, mat.entryCount());
  ASSERT_EQ(2, mat.rowCount());
  ASSERT_EQ(6, mat.colCount());
  ASSERT_EQ(5, mat.leadCol(0));
  ASSERT_EQ(0, mat.entryCountInRow(1));
  ASSERT_EQ("0: 5#101\n1:\n", mat.toString()); 
  ASSERT_TRUE(mat.emptyRow(1));

  mat.appendEntry(5, 102);
  mat.appendEntry(2000, 0); // scalar zero
  mat.rowDone(); // add a row with two entries
  ASSERT_EQ(3, mat.entryCount());
  ASSERT_EQ(3, mat.rowCount());
  ASSERT_EQ(2001, mat.colCount());
  ASSERT_EQ(5, mat.leadCol(2));
  ASSERT_EQ(2, mat.entryCountInRow(2));
  ASSERT_EQ("0: 5#101\n1:\n2: 5#102 2000#0\n", mat.toString()); 
  ASSERT_FALSE(mat.emptyRow(2));
}
