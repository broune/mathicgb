#include "stdinc.h"
#include "QuadMatrix.hpp"

#include <mathic.h>
#include <ostream>
#include <sstream>

bool QuadMatrix::debugAssertValid() const {
#ifndef MATHICGB_DEBUG
  return true;
#else
  if (!leftColumnMonomials.empty() || !rightColumnMonomials.empty()) {
    MATHICGB_ASSERT(topLeft.computeColCount() <= leftColumnMonomials.size());
    MATHICGB_ASSERT(bottomLeft.computeColCount() <=
      leftColumnMonomials.size());
    MATHICGB_ASSERT(topRight.computeColCount() <= rightColumnMonomials.size());
    MATHICGB_ASSERT(bottomRight.computeColCount() <=
    rightColumnMonomials.size());
  }

  MATHICGB_ASSERT(topLeft.rowCount() == topRight.rowCount());
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

namespace {
  class ColumnComparer {
  public:
    ColumnComparer(const PolyRing& ring): mRing(ring) {}

    typedef SparseMatrix::ColIndex ColIndex;
    typedef std::pair<monomial, ColIndex> Pair;
    bool operator()(const Pair& a, const Pair b) const {
      return mRing.monomialLT(b.first, a.first);
    }

  private:
    const PolyRing& mRing;
  };

  std::vector<SparseMatrix::ColIndex> sortColumnMonomialsAndMakePermutation(
    std::vector<monomial>& monomials,
    const PolyRing& ring
  ) {
    typedef SparseMatrix::ColIndex ColIndex;
    MATHICGB_ASSERT(monomials.size() <= std::numeric_limits<ColIndex>::max());
    const ColIndex colCount = static_cast<ColIndex>(monomials.size());
    // Monomial needs to be non-const as we are going to put these
    // monomials back into the vector of monomials which is not const.
    std::vector<std::pair<monomial, ColIndex>> columns;
    columns.reserve(colCount);
    for (ColIndex col = 0; col < colCount; ++col)
      columns.push_back(std::make_pair(monomials[col], col));
    std::sort(columns.begin(), columns.end(), ColumnComparer(ring));

    // Apply sorting permutation to monomials. This is why it is necessary to
    // copy the values in monomial out of there: in-place application of a
    // permutation is messy.
    MATHICGB_ASSERT(columns.size() == colCount);
    MATHICGB_ASSERT(monomials.size() == colCount);
    for (size_t col = 0; col < colCount; ++col) {
      MATHICGB_ASSERT(col == 0 ||
        ring.monomialLT(columns[col].first, columns[col - 1].first));
      monomials[col] = columns[col].first;
    }

    // Construct permutation of indices to match permutation of monomials
    std::vector<ColIndex> permutation(colCount);
    for (ColIndex col = 0; col < colCount; ++col) {
      // The monomial for column columns[col].second is now the
      // monomial for col, so we need the inverse map for indices.
      permutation[columns[col].second] = col;
    }

    return std::move(permutation);
  }
}

void QuadMatrix::sortColumnsLeftRightParallel(const int threadCount) {
  typedef SparseMatrix::ColIndex ColIndex;
  std::vector<ColIndex> leftPermutation;
  std::vector<ColIndex> rightPermutation;

#pragma omp parallel for num_threads(threadCount) schedule(static)
  for (OMPIndex i = 0; i < 2; ++i) {
    if (i == 0)
      leftPermutation =
        sortColumnMonomialsAndMakePermutation(leftColumnMonomials, *ring);
    else 
      rightPermutation =
        sortColumnMonomialsAndMakePermutation(rightColumnMonomials, *ring);
  }

  // todo: parallelize per block instead of per matrix.
#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (OMPIndex i = 0; i < 4; ++i) {
    if (i == 0)
      topRight.applyColumnMap(rightPermutation);
    else if (i == 1)
      bottomRight.applyColumnMap(rightPermutation);
    else if (i == 2)
      topLeft.applyColumnMap(leftPermutation);
    else {
      MATHICGB_ASSERT(i == 3);
      bottomLeft.applyColumnMap(leftPermutation);
    }
  }
}
