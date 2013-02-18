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

  QuadMatrix makeProjectionAndClear();

  const PolyRing& ring() const {return mRing;}

private:
  struct LeftRight;

  std::vector<F4ProtoMatrix*> mMatrices;

  // *** Projection of columns
  struct ColProjectTo {
    ColIndex index;
    bool isLeft;
  };
  std::vector<ColProjectTo> mColProjectTo;

  // *** Projection: Simultaneously do let/right and row permutation
  void setupRowProjection();
  struct RowProjectFrom {
    Scalar multiplyBy;
    F4ProtoMatrix::Row row;
  };
  std::vector<RowProjectFrom> mTopRowProjectFrom;
  std::vector<RowProjectFrom> mBottomRowProjectFrom;

  // *** Projection: Do left/right, then permute rows
  void setupRowProjectionLate(
    const SparseMatrix& left,
    const SparseMatrix& right
  );
  struct RowProjectFromLate {
    Scalar multiplyBy;
    RowIndex row;
  };
  void projectRows(
    SparseMatrix&& in,
    SparseMatrix& top,
    SparseMatrix& bottom
  );
  std::vector<RowProjectFromLate> mTopRows;
  std::vector<RowProjectFromLate> mBottomRows;






  std::vector<monomial> mLeftMonomials;
  std::vector<monomial> mRightMonomials;

  const PolyRing& mRing;
};

#endif
