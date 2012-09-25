#include "stdinc.h"
#include "QuadMatrixBuilder.hpp"

#include <mathic.h>
#include <sstream>

namespace {
  /// Creates a column and updates the associated data structures that
  /// are passed in. Copies mono - ownership is not taken over. The
  /// purpose of this function is to avoid code duplication. It is a
  /// template in order to avoid referring to private types of
  /// QuadMatrixBuilder.
  template<class ToMono, class ToCol>
  QuadMatrixBuilder::ColIndex createCol
  (const_monomial mono,
   SparseMatrix& top,
   SparseMatrix& bottom,
   ToMono& toMonomial,
   ToCol& toCol,
   const PolyRing& ring,
   const bool left)
  {
    MATHICGB_ASSERT(top.colCount() == bottom.colCount());
    MATHICGB_ASSERT(toMonomial.size() == bottom.colCount());
    MATHICGB_ASSERT(toCol.find(mono) == toCol.end());

    const QuadMatrixBuilder::ColIndex colCount = top.colCount();
    if (colCount == std::numeric_limits<QuadMatrixBuilder::ColIndex>::max())
      throw std::overflow_error("Too many columns in QuadMatrixBuilder");

    toMonomial.reserve(toMonomial.size() + 1);
    monomial copied = ring.allocMonomial();
    ring.monomialCopy(mono, copied);
    try {
      toCol.insert(std::make_pair
                   (copied,
                    QuadMatrixBuilder::LeftRightColIndex(colCount, left)));
    } catch (...) {
      ring.freeMonomial(copied);
      throw;
    }
    toMonomial.push_back(copied); // no throw due to reserve

    top.ensureAtLeastThisManyColumns(colCount + 1);
    bottom.ensureAtLeastThisManyColumns(colCount + 1);
    return colCount;
  }
}

QuadMatrixBuilder::ColIndex QuadMatrixBuilder::createColumnLeft
(const_monomial monomialToBeCopied) {
  return createCol
    (monomialToBeCopied,
     mTopLeft,
     mBottomLeft,
     mMonomialsLeft,
     mMonomialToCol,
     ring(),
     true);
  MATHICGB_ASSERT(mMonomialToCol.size() == leftColCount() + rightColCount());
  MATHICGB_ASSERT
    (findColumn(monomialToBeCopied).leftIndex() == leftColCount() - 1);
}

QuadMatrixBuilder::ColIndex QuadMatrixBuilder::createColumnRight
(const_monomial monomialToBeCopied) {
  return createCol
    (monomialToBeCopied,
     mTopRight,
     mBottomRight,
     mMonomialsRight,
     mMonomialToCol,
     ring(),
     false);
  MATHICGB_ASSERT(mMonomialToCol.size() == leftColCount() + rightColCount());
  MATHICGB_ASSERT
    (findColumn(monomialToBeCopied).rightIndex() == rightColCount() - 1);
}

void QuadMatrixBuilder::print(std::ostream& out) const {
  mathic::ColumnPrinter printer;
  printer.addColumn(true, "", "");
  printer.addColumn(true, " | ", "");

  // column monomials
  out << "Left columns:";
  for (ColIndex leftCol = 0; leftCol < leftColCount(); ++leftCol) {
    out << ' ';
    ring().monomialDisplay(out, monomialOfLeftCol(leftCol), false, true);
  }

  out << "\nRight columns:";
  for (ColIndex rightCol = 0; rightCol < rightColCount(); ++rightCol) {
    out << ' ';
    ring().monomialDisplay(out, monomialOfRightCol(rightCol), false, true);
  }
  out << '\n';

  // left side
  topLeft().print(printer[0]);
  printer[0] << '\n';
  bottomLeft().print(printer[0]);

  // right side
  topRight().print(printer[1]);
  printer[1] << '\n';
  bottomRight().print(printer[1]);

  out << printer;
}

// String representation intended for debugging.
std::string QuadMatrixBuilder::toString() const {
  std::ostringstream out;
  print(out);
  return out.str();
}
