#ifndef F4_MATRIX_REDUCER_GUARD
#define F4_MATRIX_REDUCER_GUARD

#include "SparseMatrix.hpp"
class QuadMatrix;
class PolyRing;

/** Class that reduces an F4 matrix represented as a QuadMatrix. The
  answer you get is the submatrix that contains new pivots. */
class F4MatrixReducer {
public:
  F4MatrixReducer(const PolyRing& ring);

  /// Reduces the lower right part of the row echelon form of matrix. Assumes
  /// that there is a permutation of the upper rows that makes the upper
  /// the upper left part of matrix upper unitriangular.
  SparseMatrix reduce(const QuadMatrix& matrix);

private:
  const SparseMatrix::Scalar mModulus;
};

#endif
