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
    if (mTopRows[lead].row != noRow && topEntryCounts[lead]<entryCount) {
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

void F4MatrixProjection::projectRowsLeftRight(
  const std::vector<RowProjectFrom>& from,
  SparseMatrix& left,
  SparseMatrix& right
) const {
  left.clear();
  right.clear();
  const auto fromEnd = from.end();
  for (auto fromIt = from.begin(); fromIt != fromEnd; ++fromIt) {
    const auto row = fromIt->row;
    MATHICGB_ASSERT(row.entryCount != 0);
    MATHICGB_ASSERT(row.scalars == 0 || row.externalScalars == 0);
    const auto modulus = static_cast<SparseMatrix::Scalar>(ring().charac());

    if (row.externalScalars != 0) {
      auto indices = row.indices;
      auto indicesEnd = row.indices + row.entryCount;
      auto scalars = row.externalScalars;
      for (; indices != indicesEnd; ++indices, ++scalars) {
        const auto scalar = static_cast<SparseMatrix::Scalar>(*scalars);
        const auto index = *indices;
        MATHICGB_ASSERT(index < mColProjectTo.size());
        const auto translated = mColProjectTo[index];
        if (translated.isLeft)
          left.appendEntry(translated.index, scalar);
        else
          right.appendEntry(translated.index, scalar);
      }
    } else {
      auto indices = row.indices;
      auto indicesEnd = row.indices + row.entryCount;
      auto scalars = row.scalars;
      for (; indices != indicesEnd; ++indices, ++scalars) {
        const auto index = *indices;
        MATHICGB_ASSERT(index < mColProjectTo.size());
        const auto translated = mColProjectTo[index];
        if (translated.isLeft)
          left.appendEntry(translated.index, *scalars);
        else
          right.appendEntry(translated.index, *scalars);
      }
    }
    const auto rowIndex = left.rowCount();
    MATHICGB_ASSERT(rowIndex == right.rowCount());
    left.rowDone();
    right.rowDone();

    if (fromIt->multiplyBy != 1) {
      MATHICGB_ASSERT(fromIt->multiplyBy != 0);
      left.multiplyRow(rowIndex, fromIt->multiplyBy, modulus);
      right.multiplyRow(rowIndex, fromIt->multiplyBy, modulus);
      MATHICGB_ASSERT(left.rowBegin(rowIndex).scalar() == 1);
    }

    MATHICGB_ASSERT(left.rowCount() == right.rowCount());
  }
}

void F4MatrixProjection::projectLeftRight(
  const std::vector<F4ProtoMatrix*>& preBlocks,
  SparseMatrix& left,
  SparseMatrix& right
) const {
  left.clear();
  right.clear();
  const auto modulus = static_cast<SparseMatrix::Scalar>(ring().charac());

  const auto end = preBlocks.end();
  for (auto it = preBlocks.begin(); it != end; ++it) {
    auto& block = **it;
    const auto rowCount = block.rowCount();
    for (SparseMatrix::RowIndex r = 0; r < rowCount; ++r) {
      const auto row = block.row(r);
      if (row.entryCount == 0)
        continue;
      MATHICGB_ASSERT(row.entryCount != 0);
      MATHICGB_ASSERT(row.scalars == 0 || row.externalScalars == 0);

      if (row.externalScalars != 0) {
        auto indices = row.indices;
        auto indicesEnd = row.indices + row.entryCount;
        auto scalars = row.externalScalars;
        for (; indices != indicesEnd; ++indices, ++scalars) {
          const auto scalar = static_cast<SparseMatrix::Scalar>(*scalars);
          const auto index = *indices;
          MATHICGB_ASSERT(index < mColProjectTo.size());
          const auto translated = mColProjectTo[index];
          if (translated.isLeft)
            left.appendEntry(translated.index, scalar);
          else
            right.appendEntry(translated.index, scalar);
        }
      } else {
        auto indices = row.indices;
        auto indicesEnd = row.indices + row.entryCount;
        auto scalars = row.scalars;
        for (; indices != indicesEnd; ++indices, ++scalars) {
          const auto index = *indices;
          MATHICGB_ASSERT(index < mColProjectTo.size());
          const auto translated = mColProjectTo[index];
          if (translated.isLeft)
            left.appendEntry(translated.index, *scalars);
          else
            right.appendEntry(translated.index, *scalars);
        }
      }
      MATHICGB_ASSERT(left.rowCount() == right.rowCount());
      left.rowDone();
      right.rowDone();
    }
  }
}

QuadMatrix F4MatrixProjection::makeProjectionAndClear() {
  QuadMatrix quadMatrix;

  if (true) {
    setupRowProjection();
    projectRowsLeftRight
      (mTopRowProjectFrom, quadMatrix.topLeft, quadMatrix.topRight);
    projectRowsLeftRight
      (mBottomRowProjectFrom, quadMatrix.bottomLeft, quadMatrix.bottomRight);
  } else {
    SparseMatrix left;
    SparseMatrix right;
    projectLeftRight(mMatrices, left, right);
    setupRowProjectionLate(left, right);
    projectRows(std::move(left), quadMatrix.topLeft, quadMatrix.bottomLeft);
    projectRows(std::move(right), quadMatrix.topRight, quadMatrix.bottomRight);
  }

  quadMatrix.ring = &ring();
  quadMatrix.leftColumnMonomials = std::move(mLeftMonomials);
  quadMatrix.rightColumnMonomials = std::move(mRightMonomials);
  return std::move(quadMatrix);
}
