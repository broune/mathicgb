#ifndef F4_MATRIX_REDUCER_GUARD
#define F4_MATRIX_REDUCER_GUARD

class QuadMatrix;
class SparseMatrix;
class PolyRing;

/** Class that reduces an F4 matrix represented as a QuadMatrix. The
  answer you get is the submatrix that contains new pivots. */
class F4MatrixReducer {
public:
  F4MatrixReducer(int threadCount);

  void reduce
  (const PolyRing& ring, QuadMatrix& matrix, SparseMatrix& newPivots);

private:
  /// this is forced to be signed to avoid warnings about signed/unsigned
  /// conversion because, perversely, MSVC 2012 does not allow unsigned
  /// for-loop indices in OpenMP.
  int mThreadCount;
};

#endif
