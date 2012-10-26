#include "stdinc.h"
#include "SparseMatrix.hpp"

#include "Poly.hpp"
#include <algorithm>

void SparseMatrix::takeRowsFrom(SparseMatrix&& matrix) {
  MATHICGB_ASSERT(matrix.colCount() <= colCount());

  if (matrix.mRows.empty())
    return;

  if (mRows.empty()) {
    const auto savedColCount = colCount();
    *this = std::move(matrix);
    mColCount = savedColCount;
    return;
  }

  Block* oldestBlock = &matrix.mBlock;
  while (oldestBlock->mPreviousBlock != 0)
    oldestBlock = oldestBlock->mPreviousBlock;

  if (mBlock.mHasNoRows) // only put mBlock in chain of blocks if non-empty
    oldestBlock->mPreviousBlock = mBlock.mPreviousBlock;
  else
    oldestBlock->mPreviousBlock = new Block(std::move(mBlock));
  mBlock = std::move(matrix.mBlock);

  mRows.insert(mRows.begin(), matrix.mRows.begin(), matrix.mRows.end());
  matrix.clear();
}

void SparseMatrix::rowToPolynomial(
  const RowIndex row,
  const std::vector<monomial>& colMonomials,
  Poly& poly
) {
  MATHICGB_ASSERT(colMonomials.size() == colCount());
  poly.setToZero();
  poly.reserve(entryCountInRow(row));
  const auto end = rowEnd(row);
  for (auto it = rowBegin(row); it != end; ++it) {
    MATHICGB_ASSERT(it.index() < colMonomials.size());
    if (it.scalar() != 0)
      poly.appendTerm(it.scalar(), colMonomials[it.index()]);
  }
  MATHICGB_ASSERT(poly.termsAreInDescendingOrder());
}

void SparseMatrix::sortRowsByIncreasingPivots() {
  SparseMatrix ordered;

  // compute pairs (pivot column index, row)
  std::vector<std::pair<SparseMatrix::ColIndex, SparseMatrix::RowIndex> > order;
  const SparseMatrix::RowIndex lRowCount = rowCount();
  const SparseMatrix::ColIndex lColCount = colCount();
  for (SparseMatrix::RowIndex row = 0; row < lRowCount; ++row) {
    if (entryCountInRow(row) == 0)
      order.push_back(std::make_pair(lColCount, row));
    else
      order.push_back(std::make_pair(rowBegin(row).index(), row));
  }

  // sort pairs by pivot column index
  std::sort(order.begin(), order.end());

  // construct ordered with pivot columns in increaing order
  ordered.clear(lColCount);
  for (size_t i = 0; i < lRowCount; ++i) {
    const auto row = order[i].second;
    const auto end = rowEnd(row);
    for (auto it = rowBegin(row); it != end; ++it)
      ordered.appendEntry(it.index(), it.scalar());
    ordered.rowDone();
  }

  *this = std::move(ordered);
}

void SparseMatrix::applyColumnMap(const std::vector<ColIndex>& colMap) {
  MATHICGB_ASSERT(colMap.size() >= colCount());
  Block* block = &mBlock;
  for (; block != 0; block = block->mPreviousBlock) {
    const auto end = block->mColIndices.end();
    for (auto it = block->mColIndices.begin(); it != end; ++it) {
      MATHICGB_ASSERT(*it < colCount());
      *it = colMap[*it];
    }
  }
}

void SparseMatrix::print(std::ostream& out) const {
  if (rowCount() == 0)
    out << "matrix with no rows\n";
  for (RowIndex row = 0; row < rowCount(); ++row) {
    out << row << ':';
    const auto end = rowEnd(row);
    for (auto it = rowBegin(row); it != end; ++it) {
      MATHICGB_ASSERT(it.index() < colCount());
      out << ' ' << it.index() << '#' << it.scalar();
    }
    out << '\n';
  }
}

std::string SparseMatrix::toString() const {
  std::ostringstream out;
  print(out);
  return out.str();
}

void SparseMatrix::appendRowAndNormalize(
  const SparseMatrix& matrix,
  const RowIndex row,
  const Scalar modulus
) {
  MATHICGB_ASSERT(row < matrix.rowCount());
  auto it = matrix.rowBegin(row);
  const auto end = matrix.rowEnd(row);
  if (it != end) {
    appendEntry(it.index(), 1);
    const Scalar lead = it.scalar();
    ++it;
    if (it != end) {
      const Scalar inverse = modularInverse(lead, modulus);
      do {
        const uint32 prod = static_cast<uint32>(inverse) * it.scalar();
        const uint16 prodMod = static_cast<uint16>(prod % modulus);
        appendEntry(it.index(), prodMod);
        ++it;
      } while (it != end);
    }
  }
  rowDone();
}

void SparseMatrix::appendRow(const SparseMatrix& matrix, const RowIndex row) {
  MATHICGB_ASSERT(row < matrix.rowCount());
  auto it = matrix.rowBegin(row);
  const auto end = matrix.rowEnd(row);
  for (; it != end; ++it)
    appendEntry(it.index(), it.scalar());
  rowDone();
}
  
SparseMatrix& SparseMatrix::operator=(const SparseMatrix& matrix) {
  clear(matrix.colCount());
  // A version that works on each block would be faster, but this is not
  // used anywhere time-critical right now. Improve this if it turns
  // up in profiling at some point.
  for (size_t row = 0; row < matrix.rowCount(); ++row)
    appendRow(matrix, row);
  return *this;
}

void SparseMatrix::swap(SparseMatrix& matrix) {
  mBlock.swap(matrix.mBlock);
  using std::swap;
  swap(mRows, matrix.mRows);
  swap(mColCount, matrix.mColCount);
}

void SparseMatrix::clear(const ColIndex newColCount) {
  Block* block = &mBlock;
  while (block != 0) {
    delete[] block->mColIndices.releaseMemory();
    delete[] block->mScalars.releaseMemory();
    Block* const tmp = block->mPreviousBlock;
    if (block != &mBlock)
      delete block;
    block = tmp;
  }
  mBlock.mPreviousBlock = 0;
  mBlock.mHasNoRows = true;
  mRows.clear();
  mColCount = newColCount;
}

void SparseMatrix::appendRowWithModulus(
  std::vector<uint64> const& v,
  const Scalar modulus
) {
  MATHICGB_ASSERT(v.size() == colCount());
  const ColIndex count = colCount();
  for (ColIndex col = 0; col < count; ++col) {
    const Scalar scalar = static_cast<Scalar>(v[col] % modulus);
    if (scalar != 0)
      appendEntry(col, scalar);
  }
  rowDone();
}

void SparseMatrix::appendRow(
  std::vector<uint64> const& v,
  const ColIndex leadCol
) {
  MATHICGB_ASSERT(v.size() == colCount());
#ifdef MATHICGB_DEBUG
  for (ColIndex col = leadCol; col < leadCol; ++col) {
    MATHICGB_ASSERT(v[col] == 0);
  }
#endif

  const ColIndex count = colCount();
  for (ColIndex col = leadCol; col < count; ++col) {
	MATHICGB_ASSERT(v[col] < std::numeric_limits<Scalar>::max());
    if (v[col] != 0)
      appendEntry(col, static_cast<Scalar>(v[col]));
  }
  rowDone();
}

void SparseMatrix::appendRowWithModulusNormalized(
  std::vector<uint64> const& v,
  const Scalar modulus
) {
  MATHICGB_ASSERT(v.size() == colCount());
  uint16 multiply = 1; 
  bool first = true;
  const ColIndex count = colCount();
  for (ColIndex col = 0; col < count; ++col) {
    Scalar scalar = static_cast<Scalar>(v[col] % modulus);
    if (scalar == 0)
      continue;
    if (first) {
      multiply = modularInverse(scalar, modulus);
      scalar = 1;
      first = false;
    } else {
      uint32 prod = static_cast<uint32>(multiply) * scalar;
      scalar = prod % modulus;
    }
    appendEntry(col, scalar);
  }
  rowDone();
}

bool SparseMatrix::appendRowWithModulusIfNonZero(
  std::vector<uint64> const& v,
  const Scalar modulus
) {
  appendRowWithModulus(v, modulus);
  MATHICGB_ASSERT(rowCount() > 0);
  if (mRows.back().empty()) {
    mRows.pop_back();
    return false;
  } else
    return true;
}

void SparseMatrix::trimLeadingZeroColumns(const ColIndex trimThisMany) {
  MATHICGB_ASSERT(trimThisMany <= colCount());
  Block* block = &mBlock;
  for (; block != 0; block = block->mPreviousBlock) {
    const auto end = block->mColIndices.end();
    for (auto it = block->mColIndices.begin(); it != end; ++it) {
      MATHICGB_ASSERT(*it >= trimThisMany);
      *it -= trimThisMany;
    }
  }
  mColCount -= trimThisMany;
}

void SparseMatrix::reserveFreeEntries(const size_t freeCount) {
  if (freeCount <= mBlock.mColIndices.capacity() - mBlock.mColIndices.size())
    return;
  // We need to copy over the pending entries, so we need space for those
  // entries on top of freeCount.
  const size_t count = freeCount + ( // todo: detect overflow for this addition
    mBlock.mHasNoRows ?
    mBlock.mColIndices.size() :
    std::distance(mRows.back().mIndicesEnd, mBlock.mColIndices.end())
  );

  auto oldBlock = new Block(std::move(mBlock));
  MATHICGB_ASSERT(mBlock.mColIndices.begin() == 0);
  MATHICGB_ASSERT(mBlock.mScalars.begin() == 0);
  MATHICGB_ASSERT(mBlock.mHasNoRows);
  MATHICGB_ASSERT(mBlock.mPreviousBlock == 0);
  mBlock.mPreviousBlock = oldBlock;

  {
    const auto begin = new ColIndex[count];
    const auto capacityEnd = begin + count;
    mBlock.mColIndices.releaseAndSetMemory(begin, begin, capacityEnd);
  }

  {
    const auto begin = new Scalar[count];
    const auto capacityEnd = begin + count;
    mBlock.mScalars.releaseAndSetMemory(begin, begin, capacityEnd);
  }

  // copy pending entries over
  if (oldBlock->mHasNoRows) {
    mBlock.mColIndices.rawAssign
      (oldBlock->mColIndices.begin(), oldBlock->mColIndices.end());
    mBlock.mScalars.rawAssign
      (oldBlock->mScalars.begin(), oldBlock->mScalars.end());
  } else {
    mBlock.mColIndices.rawAssign
      (mRows.back().mIndicesEnd, oldBlock->mColIndices.end());
    mBlock.mScalars.rawAssign
      (mRows.back().mScalarsEnd, oldBlock->mScalars.end());
  }
}

void SparseMatrix::growEntryCapacity() {
  MATHICGB_ASSERT(mBlock.mColIndices.size() == mBlock.mScalars.size());
  MATHICGB_ASSERT(mBlock.mColIndices.capacity() == mBlock.mScalars.capacity());

  // TODO: handle overflow of multiplication below
  const size_t minBlockSize = 1 << 20;
  const size_t minMultipleOfPending = 2;
  const size_t pendingCount = mBlock.mHasNoRows ?
    mBlock.mColIndices.size() :
    std::distance(mRows.back().mIndicesEnd, mBlock.mColIndices.end());
  const size_t blockSize =
    std::max(minBlockSize, pendingCount * minMultipleOfPending);

  reserveFreeEntries(blockSize);

  MATHICGB_ASSERT(mBlock.mColIndices.size() == mBlock.mScalars.size());
  MATHICGB_ASSERT(mBlock.mColIndices.capacity() == blockSize + pendingCount);
  MATHICGB_ASSERT(mBlock.mScalars.capacity() == blockSize + pendingCount);
}

std::ostream& operator<<(std::ostream& out, const SparseMatrix& matrix) {
  matrix.print(out);
  return out;
}
