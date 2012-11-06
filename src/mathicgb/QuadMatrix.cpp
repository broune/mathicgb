#include "stdinc.h"
#include "QuadMatrix.hpp"

#include <mathic.h>
#include <ostream>
#include <sstream>

bool QuadMatrix::debugAssertValid() const {
#ifndef MATHICGB_DEBUG
  return true;
#else
  MATHICGB_ASSERT(topLeft.computeColCount() <= leftColumnMonomials.size());
  MATHICGB_ASSERT(bottomLeft.computeColCount() <= leftColumnMonomials.size());
  MATHICGB_ASSERT(topLeft.rowCount() == topRight.rowCount());

  MATHICGB_ASSERT(topRight.computeColCount() <= rightColumnMonomials.size());
  MATHICGB_ASSERT(bottomRight.computeColCount() <=
    rightColumnMonomials.size());
  MATHICGB_ASSERT(bottomRight.rowCount() == bottomLeft.rowCount());
  return true;
#endif
}

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
  const auto leftColCount = leftColumnMonomials.size();
  for (ColIndex leftCol = 0; leftCol < leftColCount; ++leftCol) {
    out << ' ';
    ring->monomialDisplay(out, leftColumnMonomials[leftCol], false, true);
  }

  out << "\nRight columns:";
  const auto rightColCount = rightColumnMonomials.size();
  for (ColIndex rightCol = 0; rightCol < rightColCount; ++rightCol) {
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

size_t QuadMatrix::memoryUse() const {
  return topLeft.memoryUse() + topRight.memoryUse() +
	bottomLeft.memoryUse() + bottomRight.memoryUse();
}

size_t QuadMatrix::memoryUseTrimmed() const {
  return topLeft.memoryUseTrimmed() + topRight.memoryUseTrimmed() +
	bottomLeft.memoryUseTrimmed() + bottomRight.memoryUseTrimmed();
}

void QuadMatrix::printSizes(std::ostream& out) const {
  typedef mathic::ColumnPrinter ColPr;

  ColPr pr;
  pr.addColumn(false, " ", "");
  pr.addColumn(false, "", "");
  pr.addColumn(false, "", "");
  const char* const line = "----------";

  pr[0] << '\n';
  pr[1] << ColPr::commafy(leftColumnMonomials.size()) << "  \n";
  pr[2] << ColPr::commafy(rightColumnMonomials.size()) << "  \n";

  pr[0] << "/\n";
  pr[1] << line << "|\n";
  pr[2] << line << "\\\n";

  pr[0] << ColPr::commafy(topLeft.rowCount()) << " |\n";
  pr[1] << ColPr::commafy(topLeft.entryCount()) << " |\n";
  pr[2] << ColPr::commafy(topRight.entryCount()) << " |\n";

  pr[0] << "|\n";
  pr[1] << ColPr::bytesInUnit(topLeft.memoryUse()) << " |\n";
  pr[2] << ColPr::bytesInUnit(topRight.memoryUse()) << " |\n";

  pr[0] << "|\n";
  pr[1] << ColPr::percent(topLeft.memoryUse(), memoryUse()) << " |\n";
  pr[2] << ColPr::percent(topRight.memoryUse(), memoryUse()) << " |\n";

  pr[0] << "|\n";
  pr[1] << line << "|\n";
  pr[2] << line << "|\n";

  pr[0] << ColPr::commafy(bottomLeft.rowCount()) << " |\n";
  pr[1] << ColPr::commafy(bottomLeft.entryCount()) << " |\n";
  pr[2] << ColPr::commafy(bottomRight.entryCount()) << " |\n";

  pr[0] << "|\n";
  pr[1] << ColPr::bytesInUnit(bottomLeft.memoryUse()) << " |\n";
  pr[2] << ColPr::bytesInUnit(bottomRight.memoryUse()) << " |\n";

  pr[0] << "|\n";
  pr[1] << ColPr::percent(bottomLeft.memoryUse(), memoryUse()) << " |\n";
  pr[2] << ColPr::percent(bottomRight.memoryUse(), memoryUse()) << " |\n";

  pr[0] << "\\\n";
  pr[1] << line << "|\n";
  pr[2] << line << "/\n";

  out << '\n' << pr
	  << "Total memory: " << memoryUse() << "  ("
	  << ColPr::percent(memoryUseTrimmed(), memoryUse())
	  << " written to)\n";
}

QuadMatrix QuadMatrix::toCanonical() const {
  class RowComparer {
  public:
    RowComparer(const SparseMatrix& matrix): mMatrix(matrix) {}
    bool operator()(size_t a, size_t b) const {
      // if you need this to work for empty rows or identical leading columns
      // then update this code.
      MATHICGB_ASSERT(!mMatrix.emptyRow(a));
      MATHICGB_ASSERT(!mMatrix.emptyRow(b));
      auto itA = mMatrix.rowBegin(a);
      const auto endA = mMatrix.rowEnd(a);
      auto itB = mMatrix.rowBegin(b);
      const auto endB = mMatrix.rowEnd(b);
      for (; itA != endA; ++itA, ++itB) {
        if (itB == endB)
          return true;

        if (itA.index() > itB.index())
          return true;
        if (itA.index() < itB.index())
          return false;

        if (itA.scalar() > itB.scalar())
          return false;
        if (itA.scalar() < itB.scalar())
          return true;
      }
      return false;
    }

  private:
    const SparseMatrix& mMatrix;
  };

  const auto leftColCount = leftColumnMonomials.size();
  const auto rightColCount = rightColumnMonomials.size();

  // todo: eliminate left/right code duplication here
  QuadMatrix matrix;
  { // left side
    std::vector<size_t> rows;
    for (size_t row = 0; row < topLeft.rowCount(); ++row)
      rows.push_back(row);
    {
      RowComparer comparer(topLeft);
      std::sort(rows.begin(), rows.end(), comparer);
    }

    matrix.topLeft.clear();
    matrix.topRight.clear();
    for (size_t i = 0; i < rows.size(); ++i) {
      matrix.topLeft.appendRow(topLeft, rows[i]);
      matrix.topRight.appendRow(topRight, rows[i]);
    }
  }
  { // right side
    std::vector<size_t> rows;
    for (size_t row = 0; row < bottomLeft.rowCount(); ++row)
      rows.push_back(row);
    {
      RowComparer comparer(bottomLeft);
      std::sort(rows.begin(), rows.end(), comparer);
    }

    matrix.bottomLeft.clear();
    matrix.bottomRight.clear();
    for (size_t i = 0; i < rows.size(); ++i) {
      matrix.bottomLeft.appendRow(bottomLeft, rows[i]);
      matrix.bottomRight.appendRow(bottomRight, rows[i]);
    }
  }

  matrix.leftColumnMonomials = leftColumnMonomials;
  matrix.rightColumnMonomials = rightColumnMonomials;
  matrix.ring = ring;
  
  return std::move(matrix);
}

std::ostream& operator<<(std::ostream& out, const QuadMatrix& qm) {
  qm.print(out);
  return out;
}
