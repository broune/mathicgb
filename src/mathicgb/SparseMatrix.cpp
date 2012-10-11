#include "stdinc.h"
#include "SparseMatrix.hpp"

#include "Poly.hpp"
#include <algorithm>

std::ostream& operator<<(std::ostream& out, const SparseMatrix& matrix) {
  matrix.print(out);
  return out;
}

void SparseMatrix::rowToPolynomial
(RowIndex row, std::vector<monomial> colMonomials, Poly& poly) {
  MATHICGB_ASSERT(colMonomials.size() == colCount());
  poly.setToZero();
  auto end = rowEnd(row);
  poly.reserve(entryCountInRow(row));
  for (auto it = rowBegin(row); it != end; ++it) {
    MATHICGB_ASSERT(it.index() < colMonomials.size());
    if (it.scalar() != 0)
      poly.appendTerm(it.scalar(), colMonomials[it.index()]);
  }
}

void SparseMatrix::sortRowsByIncreasingPivots() {
  SparseMatrix ordered;

  // compute pairs (pivot column index, row)
  std::vector<std::pair<SparseMatrix::ColIndex, SparseMatrix::RowIndex> > order;
  const SparseMatrix::RowIndex lRowCount = rowCount();
  const SparseMatrix::ColIndex lColCount = colCount();
  for (SparseMatrix::RowIndex row = 0; row < lRowCount; ++row) {
    if (entryCountInRow(row) == 0)
      order.push_back(std::make_pair(lColCount, row));
    else
      order.push_back(std::make_pair(rowBegin(row).index(), row));
  }

  // sort pairs by pivot column index
  std::sort(order.begin(), order.end());

  // construct ordered with pivot columns in increaing order
  ordered.clear(lColCount);
  for (size_t i = 0; i < lRowCount; ++i) {
    const SparseMatrix::RowIndex row = order[i].second;
    SparseMatrix::RowIterator it = rowBegin(row);
    SparseMatrix::RowIterator end = rowEnd(row);
    for (; it != end; ++it)
      ordered.appendEntry(it.index(), it.scalar());
    ordered.rowDone();
  }

  *this = std::move(ordered);
}
