// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_QUAD_MATRIX_GUARD
#define MATHICGB_QUAD_MATRIX_GUARD

#include "PolyRing.hpp"
#include "SparseMatrix.hpp"
#include <vector>
#include <string>
#include <ostream>

MATHICGB_NAMESPACE_BEGIN

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

  void printStatistics(std::ostream& out) const;

  size_t memoryUse() const;
  size_t memoryUseTrimmed() const;

  /// Shows whole matrix in a string. Useful for debugging.
  std::string toString() const;

  /// Return the combined number of non-zero entries.
  size_t entryCount() const;

  /// Return the combined number of left and right columns.
  size_t rowCount() const;

  /// Return the number of left columns.
  SparseMatrix::ColIndex computeLeftColCount() const;

  /// Return the number of right columns.
  SparseMatrix::ColIndex computeRightColCount() const;

  void write(SparseMatrix::Scalar modulus, FILE* file) const;

  /// Read a matrix from file into *this. Return the modulus from file.
  /// This method clears the column monomials and the ring pointer.
  SparseMatrix::Scalar read(FILE* file);

  /// Sort the left columns to be in decreasing order according to the monomial
  /// order from the ring. The operation is done in parallel.
  void sortColumnsLeftRightParallel();

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

MATHICGB_NAMESPACE_END

#endif
