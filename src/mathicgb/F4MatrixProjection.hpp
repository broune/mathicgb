#ifndef MATHICGB_F4_MATRIX_PROJECTION_GUARD
#define MATHICGB_F4_MATRIX_PROJECTION_GUARD

#include "QuadMatrix.hpp"
#include "SparseMatrix.hpp"
#include "F4ProtoMatrix.hpp"
#include "MonomialMap.hpp"
#include <vector>

class F4MatrixProjection {
public:
  typedef SparseMatrix::ColIndex ColIndex;

  F4MatrixProjection() {}

  QuadMatrix project(
    std::vector<char>&& isColumnToLeft,
    MonomialMap<ColIndex>&& map,
    std::vector<F4ProtoMatrix*>&& matrix
  );

private:
};

#endif
