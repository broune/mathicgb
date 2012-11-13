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
    topLeft(std::move(matrix.topLeft)),
    topRight(std::move(matrix.topRight)),
    bottomLeft(std::move(matrix.bottomLeft)),
    bottomRight(std::move(matrix.bottomRight)),
    leftColumnMonomials(std::move(matrix.leftColumnMonomials)),
    rightColumnMonomials(std::move(matrix.rightColumnMonomials)),
    ring(matrix.ring)
  {}
  QuadMatrix& operator=(QuadMatrix&& matrix) {
    this->~QuadMatrix();
    new (this) QuadMatrix(std::move(matrix));
    return *this;
  }

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

  /// Prints the sizes of the matrix out, in terms of number of rows and
  /// columns and how many non-zero entries in each submatrix.
  void printSizes(std::ostream& out) const;

  size_t memoryUse() const;
  size_t memoryUseTrimmed() const;

  /// Shows whole matrix in a string. Useful for debugging.
  std::string toString() const;

  /// Sort the left columns to be in decreasing order according to the monomial
  /// order from the ring. The operation is parallel using up to threadCount
  /// threads.
  void sortColumnsLeftRightParallel(int threadCount);

  /// Makes a copy of this matrix whose rows are sorted in some canonical way.
  /// TODO: Actually only coarsely sorts the top rows right now.
  QuadMatrix toCanonical() const;

  /// Asserts internal invariants if asserts are turned on.
  bool debugAssertValid() const;

private:
  QuadMatrix(const QuadMatrix&); // not available
  void operator=(const QuadMatrix&); // not available
};

std::ostream& operator<<(std::ostream& out, const QuadMatrix& qm);

#endif
