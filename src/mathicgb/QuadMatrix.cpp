#include "stdinc.h"
#include "QuadMatrix.hpp"

#include <mathic.h>
#include <tbb/tbb.h>
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

size_t QuadMatrix::rowCount() const {
  return topLeft.rowCount() + bottomLeft.rowCount();
}

SparseMatrix::ColIndex QuadMatrix::computeLeftColCount() const {
  if (!leftColumnMonomials.empty()) {
    MATHICGB_ASSERT(
      leftColumnMonomials.size() <=
      std::numeric_limits<SparseMatrix::ColIndex>::max()
    );
    return static_cast<SparseMatrix::ColIndex>(leftColumnMonomials.size());
  }
  return std::max(topLeft.computeColCount(), bottomLeft.computeColCount());
}

SparseMatrix::ColIndex QuadMatrix::computeRightColCount() const {
  if (!rightColumnMonomials.empty()) {
    MATHICGB_ASSERT(
      rightColumnMonomials.size() <=
      std::numeric_limits<SparseMatrix::ColIndex>::max()
    );
    return static_cast<SparseMatrix::ColIndex>(rightColumnMonomials.size());
  }
  return std::max(topRight.computeColCount(), bottomRight.computeColCount());
}

size_t QuadMatrix::entryCount() const {
  return
    topLeft.entryCount() + topRight.entryCount() +
    bottomLeft.entryCount() + bottomRight.entryCount();
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

void QuadMatrix::printStatistics(std::ostream& out) const {
  typedef mathic::ColumnPrinter ColPr;

  ColPr pr;
  pr.addColumn(false, " ", "");
  pr.addColumn(false, "", "");
  pr.addColumn(false, "", "");
  pr.addColumn(false, "", "");
  pr.addColumn(true, "", "");

  const auto totalMemory = memoryUse();

  auto printDataCol = [&](
    std::ostream& out,
    const SparseMatrix& top,
    const SparseMatrix& bottom,
    const SparseMatrix::ColIndex colCount
  ) {
    auto printDataCell = [&](const SparseMatrix& matrix) {
      const auto entryCount = matrix.entryCount();
      const uint64 area =
        static_cast<uint64>(matrix.rowCount()) * static_cast<uint64>(colCount);
      const auto memory = matrix.memoryUse();

      out << ColPr::withSIPrefix(entryCount) << " -"
        << ColPr::percentFixed(entryCount, area) << " \n"
        << ColPr::bytesInUnit(memory) << " -"
        << ColPr::percentFixed(matrix.memoryUseTrimmed(), memory) << " \n"
        << ColPr::percent(memory, totalMemory) << " \n";
    };

    out << ColPr::commafy(colCount) << " \n";
    const char* const line = "------------------\n";
    out << line;
    printDataCell(top);
    out << line;
    printDataCell(bottom);
    out << line;
  };

  pr[0] << "\n/\n" << ColPr::commafy(topLeft.rowCount())
    << " |\nrows |\n|\n|\n"
    << ColPr::commafy(bottomLeft.rowCount())
    << " |\nrows |\n|\n\\\n";
  printDataCol(pr[1], topLeft, bottomLeft, computeLeftColCount());
  pr[2] << " \n|\n|\n|\n|\n|\n|\n|\n|\n|\n";
  printDataCol(pr[3], topRight, bottomRight, computeRightColCount());

  const char* const legend =
    "| non-zero (density)\n| memory (used)\n| of total memory\n";
  pr[4] << "  columns\n\\\n" << legend << "|\n" << legend << "/\n";

  out << '\n' << pr
    << "       Total memory: " << ColPr::bytesInUnit(memoryUse()) << " ("
	<< ColPr::percent(memoryUseTrimmed(), totalMemory)
	<< " used)\n\n";
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

void QuadMatrix::sortColumnsLeftRightParallel() {
  typedef SparseMatrix::ColIndex ColIndex;
  std::vector<ColIndex> leftPermutation;
  std::vector<ColIndex> rightPermutation;
  
  tbb::parallel_for(0, 2, 1, [&](int i) {
    if (i == 0)
      leftPermutation =
        sortColumnMonomialsAndMakePermutation(leftColumnMonomials, *ring);
    else 
      rightPermutation =
        sortColumnMonomialsAndMakePermutation(rightColumnMonomials, *ring);
  });

  tbb::parallel_for(0, 4, 1, [&](int i) {
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
  });
}

void QuadMatrix::write(
  const SparseMatrix::Scalar modulus,
  FILE* file
) const {
  MATHICGB_ASSERT(file != 0);
  topLeft.write(modulus, file);
  topRight.write(modulus, file);
  bottomLeft.write(modulus, file);
  bottomRight.write(modulus, file);
}

SparseMatrix::Scalar QuadMatrix::read(FILE* file) {
  MATHICGB_ASSERT(file != 0);

  leftColumnMonomials.clear();
  rightColumnMonomials.clear();
  ring = 0;

  const auto topLeftModulus = topLeft.read(file);
  const auto topRightModulus = topRight.read(file);
  const auto bottomLeftModulus = bottomLeft.read(file);
  const auto bottomRightModulus = bottomRight.read(file);

  // todo: this should throw some kind of invalid format exception instead of
  // these asserts.
  MATHICGB_ASSERT(topLeftModulus == topRightModulus);
  MATHICGB_ASSERT(bottomLeftModulus == topRightModulus);
  MATHICGB_ASSERT(bottomRightModulus == topRightModulus);
  MATHICGB_ASSERT(debugAssertValid());

  return topLeftModulus;
}
