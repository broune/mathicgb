#include "stdinc.h"
#include "SparseMatrix.hpp"

#include "Poly.hpp"

std::ostream& operator<<(std::ostream& out, const SparseMatrix& matrix) {
  matrix.print(out);
  return out;
}

void SparseMatrix::rowToPolynomial
(RowIndex row, std::vector<monomial> colMonomials, Poly& poly) {
  MATHICGB_ASSERT(colMonomials.size() == colCount());
  poly.setToZero();
  auto end = rowEnd(row);
  for (auto it = rowBegin(row); it != end; ++it) {
    MATHICGB_ASSERT(it.index() < colMonomials.size());
    if (it.scalar() != 0)
      poly.appendTerm(it.scalar(), colMonomials[it.index()]);
  }
}
