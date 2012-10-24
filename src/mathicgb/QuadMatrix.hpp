#ifndef MATHICGB_QUAD_MATRIX_GUARD
#define MATHICGB_QUAD_MATRIX_GUARD

#include "PolyRing.hpp"
#include "SparseMatrix.hpp"
#include <vector>
#include <string>
#include <ostream>
class ostream;

/** Represents a matrix composed of 4 sub-matrices that fit together
  into one matrix divided into top left, top right, bottom left and
  bottom right. This is a convenient representation of the matrices
  encountered in the F4 polynomial reduction algorithm.

  So far this is just a data class used as the output of a
  QuadMatrixBuilder or F4MatrixBuilder. */
class QuadMatrix {
public:
  QuadMatrix() {}
  QuadMatrix(QuadMatrix&& matrix):
    topLeft(matrix.topLeft),
    topRight(matrix.topRight),
    bottomLeft(matrix.bottomLeft),
    bottomRight(matrix.bottomRight),
    leftColumnMonomials(matrix.leftColumnMonomials),
    rightColumnMonomials(matrix.rightColumnMonomials),
    ring(matrix.ring)
  {}

  SparseMatrix topLeft; 
  SparseMatrix topRight;
  SparseMatrix bottomLeft;
  SparseMatrix bottomRight;
  std::vector<monomial> leftColumnMonomials;
  std::vector<monomial> rightColumnMonomials;
  const PolyRing* ring;

  /// Prints whole matrix to out in human-readable format. Useful for
  /// debugging.
  void print(std::ostream& out) const;

  /// Shows whole matrix in a string. Useful for debugging.
  std::string toString() const;

  /// Makes a copy of this matrix whose rows are sorted in some canonical way.
  /// TODO: Actually only coarsely sorts the top rows right now.
  QuadMatrix toCanonical() const;

#ifdef MATHICGB_DEBUG
  bool debugAssertValid() const;
#endif

private:
  QuadMatrix(const QuadMatrix&); // not available
  void operator=(const QuadMatrix&); // not available
};

std::ostream& operator<<(std::ostream& out, const QuadMatrix& qm);

#endif
