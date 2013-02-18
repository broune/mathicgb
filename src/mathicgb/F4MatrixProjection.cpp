#include "stdinc.h"
#include "F4MatrixProjection.hpp"

#include "ScopeExit.hpp"

F4MatrixProjection::F4MatrixProjection(
  const PolyRing& ring,
  ColIndex colCount
):
  mRing(ring),
  mColProjectTo(colCount)
{}

void F4MatrixProjection::addColumn(
  const ColIndex projectFrom,
  const const_monomial mono,
  const bool isLeft
) {
  MATHICGB_ASSERT(projectFrom < mColProjectTo.size());
  MATHICGB_ASSERT
    (mLeftMonomials.size() + mRightMonomials.size() < mColProjectTo.size());

  auto monoCopy = mRing.allocMonomial();
  MATHICGB_SCOPE_EXIT(monoGuard) {mRing.freeMonomial(monoCopy);};
  mRing.monomialCopy(mono, monoCopy);

  auto& projected = mColProjectTo[projectFrom];
  if (isLeft) {
    projected.isLeft = true;
    projected.index = static_cast<ColIndex>(mLeftMonomials.size());
    mLeftMonomials.push_back(monoCopy);
  } else {
    projected.isLeft = false;
    projected.index = static_cast<ColIndex>(mRightMonomials.size());
    mRightMonomials.push_back(monoCopy);
  }

  monoGuard.dismiss();
}

void F4MatrixProjection::setupRowProjection() {
  mTopRowProjectFrom.resize(mLeftMonomials.size());

  const auto noReducer = std::numeric_limits<RowIndex>::max();
  F4ProtoMatrix::Row noRow = {};
  noRow.indices = 0;
  const auto noCol = std::numeric_limits<SparseMatrix::ColIndex>::max();
  const auto modulus = static_cast<SparseMatrix::Scalar>(ring().charac());

  const auto end = mMatrices.end();
  for (auto it = mMatrices.begin(); it != end; ++it) {
    auto& matrix = **it;
    const auto rowCount = matrix.rowCount();
    for (RowIndex r = 0; r < rowCount; ++r) {
      const auto row = matrix.row(r);
      if (row.entryCount == 0)
        continue;

      // Determine leading (minimum index) left entry.
      const auto lead = [&] {
        for (SparseMatrix::ColIndex col = 0; col < row.entryCount; ++col) {
          auto const translated = mColProjectTo[row.indices[col]];
          if (translated.isLeft)
            return std::make_pair(col, translated.index);
        }
        return std::make_pair(noCol, noCol);
      }();
      // Decide if this should be a reducer or reducee row.
      if (lead.second == noCol) {
        const RowProjectFrom p = {1, row};
        mBottomRowProjectFrom.push_back(p); // no left entries
        continue;
      }
      MATHICGB_ASSERT(row.scalars != 0 || row.externalScalars != 0);

      const auto reducer = mTopRowProjectFrom[lead.second].row;
      if (
        reducer.entryCount != 0 && // already have a reducer and...
        reducer.entryCount < row.entryCount // ...it is sparser/better
      ) {
        const RowProjectFrom p = {1, row};
        mBottomRowProjectFrom.push_back(p);
      } else {
        if (reducer.entryCount != 0) {
          const RowProjectFrom p = {1, reducer};
          mBottomRowProjectFrom.push_back(p);
        }
        const auto leadScalar = row.scalars != 0 ? row.scalars[lead.first] :
          static_cast<SparseMatrix::Scalar>
            (row.externalScalars[lead.first]);
        MATHICGB_ASSERT(leadScalar != 0);
        const auto inverse = leadScalar == 1 ?
          1 : modularInverse(leadScalar, modulus);
        const RowProjectFrom p = {inverse, row};
        mTopRowProjectFrom[lead.second] = p;
      }
    }
  }
#ifdef MATHICGB_DEBUG
  for (size_t  i = 0; i < mTopRowProjectFrom.size(); ++i) {
    const auto p = mTopRowProjectFrom[i];
    MATHICGB_ASSERT(p.row.entryCount > 0);
    MATHICGB_ASSERT(p.multiplyBy != 0);

    // Find leading left entry.
    ColIndex col = 0;
    while (!mColProjectTo[p.row.indices[col]].isLeft) {
      ++col;
      MATHICGB_ASSERT(col < p.row.entryCount);
    }
    const auto projected = mColProjectTo[p.row.indices[col]];

    // The leading left entry of row i should be in column i.
    MATHICGB_ASSERT(projected.index == i);

    // After multiplication, the leading left scalar should be 1.
    const auto leadScalar = p.row.scalars != 0 ?
      p.row.scalars[col] : static_cast<Scalar>(p.row.externalScalars[col]);
    MATHICGB_ASSERT(modularProduct(leadScalar, p.multiplyBy, modulus) == 1);
  }
  for (size_t  i = 0; i < mBottomRowProjectFrom.size(); ++i) {
    const auto p = mBottomRowProjectFrom[i];
    MATHICGB_ASSERT(p.row.entryCount > 0);
    MATHICGB_ASSERT(p.multiplyBy != 0);
  }
#endif
}

void F4MatrixProjection::projectRows(
  SparseMatrix&& in,
  SparseMatrix& top,
  SparseMatrix& bottom
) {
  const auto modulus = static_cast<Scalar>(ring().charac());

  top.clear();
  bottom.clear();

  const auto rowCountTop =
    static_cast<SparseMatrix::RowIndex>(mTopRows.size());
  for (SparseMatrix::RowIndex toRow = 0; toRow < rowCountTop; ++toRow) {
    top.appendRow(in, mTopRows[toRow].row);
    if (mTopRows[toRow].multiplyBy != 1)
      top.multiplyRow(toRow, mTopRows[toRow].multiplyBy, modulus);
  }

  const auto rowCountBottom =
    static_cast<SparseMatrix::RowIndex>(mBottomRows.size());
  for (SparseMatrix::RowIndex toRow = 0; toRow < rowCountBottom; ++toRow) {
    bottom.appendRow(in, mBottomRows[toRow].row);
    if (mBottomRows[toRow].multiplyBy != 1)
        bottom.multiplyRow(toRow, mBottomRows[toRow].multiplyBy, modulus);
  }

  in.clear();
}

void F4MatrixProjection::setupRowProjectionLate(
  const SparseMatrix& left,
  const SparseMatrix& right
) {
  const auto leftColCount = static_cast<ColIndex>(mLeftMonomials.size());
  const auto modulus = static_cast<Scalar>(ring().charac());
  const auto noRow = std::numeric_limits<ColIndex>::max();
  {
    const RowProjectFromLate init = {0, noRow};
    mTopRows.resize(leftColCount, init);
  }

  MATHICGB_ASSERT(left.computeColCount() == leftColCount);
  MATHICGB_ASSERT(left.rowCount() >= leftColCount);
  MATHICGB_ASSERT(left.rowCount() == right.rowCount());

  std::vector<SparseMatrix::ColIndex> topEntryCounts(leftColCount);

  const auto rowCount = left.rowCount();
  for (SparseMatrix::RowIndex row = 0; row < rowCount; ++row) {
    const auto leftEntryCount = left.entryCountInRow(row);
    const auto entryCount = leftEntryCount + right.entryCountInRow(row);
    MATHICGB_ASSERT(entryCount >= leftEntryCount); // no overflow
    if (entryCount == 0)
      continue; // ignore zero rows
    if (leftEntryCount == 0) {
      RowProjectFromLate p = {1, row};
      mBottomRows.push_back(p); // can't be top/reducer
      continue;
    }
    const auto lead = left.rowBegin(row).index();
    if (mTopRows[lead].row != noRow && topEntryCounts[lead] < entryCount) {
      RowProjectFromLate p = {1, row};
      mBottomRows.push_back(p); // other reducer better
    } else {
      if (mTopRows[lead].row != noRow) {
        RowProjectFromLate p = {1, mTopRows[lead].row};
        mBottomRows.push_back(p);
      }
      topEntryCounts[lead] = entryCount;
      mTopRows[lead].row = row; // do scalar .first later
    }
  }

  for (SparseMatrix::RowIndex r = 0; r < leftColCount; ++r) {
    const auto row = mTopRows[r].row;
    MATHICGB_ASSERT(row != noRow);
    MATHICGB_ASSERT(left.entryCountInRow(row) > 0);
    MATHICGB_ASSERT(left.rowBegin(row).index() == r);
    MATHICGB_ASSERT(left.rowBegin(row).scalar() != 0);
    MATHICGB_ASSERT(topEntryCounts[r] ==
      left.entryCountInRow(row) + right.entryCountInRow(row));

    const auto leadScalar = left.rowBegin(row).scalar();
    mTopRows[r].multiplyBy = leadScalar == 1 ? 1 : // 1 is the common case
      modularInverse(leadScalar, modulus);
  }

#ifdef MATHICGB_DEBUG
  for (SparseMatrix::RowIndex r = 0; r < mBottomRows.size(); ++r) {
    const auto row = mBottomRows[r].row;
    MATHICGB_ASSERT(
      left.entryCountInRow(row) + right.entryCountInRow(row) > 0);
    MATHICGB_ASSERT(mBottomRows[r].multiplyBy == 1);
  }
#endif
}

// Utility class for building a left/right projection.
struct F4MatrixProjection::LeftRight {
  typedef F4ProtoMatrix::ExternalScalar ExternalScalar;
  typedef F4ProtoMatrix::Row Row;

  LeftRight(
    const std::vector<ColProjectTo>& colProjectTo,
    const PolyRing& ring,
    const size_t quantum
  ):
    mColProjectTo(colProjectTo),
    mModulus(static_cast<Scalar>(ring.charac())),
    mLeft(quantum),
    mRight(quantum)
  {
    MATHICGB_ASSERT(ring.charac() < std::numeric_limits<Scalar>::max());
    mLeft.clear();
    mRight.clear();
  }

  void appendRowsPermuted(const std::vector<RowProjectFrom>& rows) {
    const auto end = rows.end();
    for (auto it = rows.begin(); it != end; ++it)
      appendRow(it->row, it->multiplyBy);
  }

  void appendRows(const std::vector<F4ProtoMatrix*>& preBlocks) {
    const auto end = preBlocks.end();
    for (auto it = preBlocks.begin(); it != end; ++it) {
      auto& block = **it;
      const auto rowCount = block.rowCount();
      for (SparseMatrix::RowIndex r = 0; r < rowCount; ++r) {
        const auto row = block.row(r);
        if (row.entryCount > 0)
          appendRow(row);
      }
    }
  }

  void appendRow(const Row& row) {
    MATHICGB_ASSERT(row.entryCount > 0); // could be OK, but unexpected
    MATHICGB_ASSERT(row.scalars == 0 || row.externalScalars == 0);

    const auto indicesEnd = row.indices + row.entryCount;
    if (row.externalScalars != 0)
      appendRow(row.indices, indicesEnd, row.externalScalars);
    else
      appendRow(row.indices, indicesEnd, row.scalars);
  }

  void appendRow(const Row& row, Scalar multiplyBy) {
    MATHICGB_ASSERT(multiplyBy != 0);
    appendRow(row);
    if (multiplyBy != 1) {
      const auto rowIndex = mLeft.rowCount() - 1;
      mLeft.multiplyRow(rowIndex, multiplyBy, mModulus);
      mRight.multiplyRow(rowIndex, multiplyBy, mModulus);
    }
  }

  template<class IndexIter, class ScalarIter>
  void appendRow(
    IndexIter indices,
    const IndexIter indicesEnd,
    ScalarIter scalars
  ) {
    for (; indices != indicesEnd; ++indices, ++scalars)
      appendEntry(*indices, *scalars);
    rowDone();
  }

  void appendEntry(const ColIndex projectMe, const Scalar scalar) {
    MATHICGB_ASSERT(scalar < mModulus);
    MATHICGB_ASSERT(mLeft.rowCount() == mRight.rowCount());
    MATHICGB_ASSERT(projectMe < mColProjectTo.size());
    const auto projected = mColProjectTo[projectMe];
    if (projected.isLeft)
      mLeft.appendEntry(projected.index, scalar);
    else
      mRight.appendEntry(projected.index, scalar);
  }

  void appendEntry(const ColIndex projectMe, const ExternalScalar scalar) {
    MATHICGB_ASSERT(scalar <= std::numeric_limits<Scalar>::max());
    appendEntry(projectMe, static_cast<Scalar>(scalar));
  }

  void rowDone() {
    MATHICGB_ASSERT(mLeft.rowCount() == mRight.rowCount());
    mLeft.rowDone();
    mRight.rowDone();
  };

  const SparseMatrix& left() const {return mLeft;}
  const SparseMatrix& right() const {return mRight;}

  SparseMatrix moveLeft() {return std::move(mLeft);}
  SparseMatrix moveRight() {return std::move(mRight);}

private:
  const std::vector<ColProjectTo>& mColProjectTo;
  const Scalar mModulus;

  SparseMatrix mLeft;
  SparseMatrix mRight;
};

QuadMatrix F4MatrixProjection::makeProjectionAndClear() {
  QuadMatrix qm;

  const size_t quantum = 0; // todo: set quantum from parameter

  if (true) {
    setupRowProjection();
    {
      LeftRight top(mColProjectTo, ring(), 0);
      top.appendRowsPermuted(mTopRowProjectFrom);
      qm.topLeft = top.moveLeft();
      qm.topRight = top.moveRight();
    }
    {
      LeftRight bottom(mColProjectTo, ring(), 0);
      bottom.appendRowsPermuted(mBottomRowProjectFrom);
      qm.bottomLeft = bottom.moveLeft();
      qm.bottomRight = bottom.moveRight();
    }
  } else {
    LeftRight lr(mColProjectTo, ring(), 0);
    lr.appendRows(mMatrices);
    setupRowProjectionLate(lr.left(), lr.right());
    projectRows(lr.moveLeft(), qm.topLeft, qm.bottomLeft);
    projectRows(lr.moveRight(), qm.topRight, qm.bottomRight);
  }

  qm.ring = &ring();
  qm.leftColumnMonomials = std::move(mLeftMonomials);
  qm.rightColumnMonomials = std::move(mRightMonomials);
  return std::move(qm);
}
