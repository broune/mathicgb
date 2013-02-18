#ifndef MATHICGB_F4_MATRIX_PROJECTION_GUARD
#define MATHICGB_F4_MATRIX_PROJECTION_GUARD

#include "QuadMatrix.hpp"
#include "SparseMatrix.hpp"
#include "F4ProtoMatrix.hpp"
#include "MonomialMap.hpp"
#include "ScopeExit.hpp"
#include <vector>

class F4MatrixProjection {
public:
  typedef SparseMatrix::RowIndex RowIndex;
  typedef SparseMatrix::ColIndex ColIndex;
  typedef SparseMatrix::Scalar Scalar;

  F4MatrixProjection(const PolyRing& ring, ColIndex colCount);

  void addProtoMatrix(F4ProtoMatrix&& matrix) {mMatrices.push_back(&matrix);}

  // No reference to mono is retained.
  void addColumn(ColIndex index, const_monomial mono, const bool isLeft);

  QuadMatrix makeAndClear(const size_t quantum);

  const PolyRing& ring() const {return mRing;}

private:
  QuadMatrix makeAndClearOneStep(const size_t quantum);
  QuadMatrix makeAndClearTwoStep(const size_t quantum);

  // Utility class for building a left/right projection.
  class LeftRight;

  // Utility class for building a top/bottom projection.
  template<class Row>
  class TopBottom;

  // This is for projection of columns
  struct ColProjectTo {
    ColIndex index;
    bool isLeft;
  };
  std::vector<ColProjectTo> mColProjectTo;

  std::vector<F4ProtoMatrix*> mMatrices;
  std::vector<monomial> mLeftMonomials;
  std::vector<monomial> mRightMonomials;
  const PolyRing& mRing;
};

#endif
