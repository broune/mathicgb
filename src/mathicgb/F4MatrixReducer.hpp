#ifndef F4_MATRIX_REDUCER_GUARD
#define F4_MATRIX_REDUCER_GUARD

#include "SparseMatrix.hpp"
class QuadMatrix;
class PolyRing;

/** Class that reduces an F4 matrix represented as a QuadMatrix. The
  answer you get is the submatrix that contains new pivots. */
class F4MatrixReducer {
public:
  F4MatrixReducer(const PolyRing& ring, int threadCount);

  SparseMatrix reduce(QuadMatrix& matrix);

private:
  /// this is forced to be signed to avoid warnings about signed/unsigned
  /// conversion because, perversely, MSVC 2012 does not allow unsigned
  /// for-loop indices in OpenMP.
  const SparseMatrix::Scalar mModulus;
  const int mThreadCount;
};

#endif
