#include "stdinc.h"
#include "QuadMatrix.hpp"

#include <ostream>
#include <sstream>
#include <mathic.h>

#ifdef MATHICGB_DEBUG
bool QuadMatrix::debugAssertValid() const {
  MATHICGB_ASSERT(topLeft.colCount() == bottomLeft.colCount());
  MATHICGB_ASSERT(topLeft.rowCount() == topRight.rowCount());
  MATHICGB_ASSERT(topLeft.colCount() == leftColumnMonomials.size());

  MATHICGB_ASSERT(bottomRight.colCount() == topRight.colCount());
  MATHICGB_ASSERT(bottomRight.rowCount() == bottomLeft.rowCount());
  MATHICGB_ASSERT(bottomRight.colCount() == rightColumnMonomials.size());   
  return true;
}
#endif

void QuadMatrix::print(std::ostream& out) const {
  MATHICGB_ASSERT(debugAssertValid());

  // @todo: fix the code duplication here from QuadMatrixBuilder's
  // string code.

  typedef SparseMatrix::ColIndex ColIndex;
  mathic::ColumnPrinter printer;
  printer.addColumn(true, "", "");
  printer.addColumn(true, " | ", "");

  // column monomials
  out << "Left columns:";
  for (ColIndex leftCol = 0; leftCol < topLeft.colCount(); ++leftCol) {
    out << ' ';
    ring->monomialDisplay(out, leftColumnMonomials[leftCol], false, true);
  }

  out << "\nRight columns:";
  for (ColIndex rightCol = 0; rightCol < topRight.colCount(); ++rightCol) {
    out << ' ';
    ring->monomialDisplay(out, rightColumnMonomials[rightCol], false, true);
  }
  out << '\n';

  // left side
  topLeft.print(printer[0]);
  printer[0] << '\n';
  bottomLeft.print(printer[0]);

  // right side
  topRight.print(printer[1]);
  printer[1] << '\n';
  bottomRight.print(printer[1]);

  out << printer;
}

std::string QuadMatrix::toString() const {
  std::ostringstream out;
  print(out);
  return out.str();
}

QuadMatrix QuadMatrix::toCanonical() const {
  std::vector<size_t> rows;
  for (size_t row = 0; row < topLeft.rowCount(); ++row)
    rows.push_back(row);
  class RowComparer {
  public:
    RowComparer(const SparseMatrix& matrix): mMatrix(matrix) {}
    bool operator()(size_t a, size_t b) const {
      // if you need this to work for empty rows or identical leading columns
      // then update this code.
      MATHICGB_ASSERT(!mMatrix.emptyRow(a));
      MATHICGB_ASSERT(!mMatrix.emptyRow(b));
      return mMatrix.leadCol(a) > mMatrix.leadCol(b);
    }

  private:
    const SparseMatrix& mMatrix;
  };
  {
    RowComparer comparer(topLeft);
    std::sort(rows.begin(), rows.end(), comparer);
  }

  QuadMatrix matrix;
  matrix.topLeft.clear(topLeft.colCount());
  matrix.topRight.clear(topRight.colCount());
  for (size_t i = 0; i < rows.size(); ++i) {
    matrix.topLeft.appendRow(topLeft, rows[i]);
    matrix.topRight.appendRow(topRight, rows[i]);
  }

  matrix.bottomLeft = bottomLeft;
  matrix.bottomRight = bottomRight;
  matrix.leftColumnMonomials = leftColumnMonomials;
  matrix.rightColumnMonomials = rightColumnMonomials;
  matrix.ring = ring;
  
  return std::move(matrix);
}

std::ostream& operator<<(std::ostream& out, const QuadMatrix& qm) {
  qm.print(out);
  return out;
}
