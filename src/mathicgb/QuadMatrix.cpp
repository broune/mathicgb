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

// String representation intended for debugging.
std::string QuadMatrix::toString() const {
  std::ostringstream out;
  print(out);
  return out.str();
}
