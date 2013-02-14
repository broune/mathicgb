#include "stdinc.h"
#include "F4MatrixProjection.hpp"

namespace {
  class LeftRightProjection {
  public:
    typedef SparseMatrix::ColIndex ColIndex;

    LeftRightProjection(
      const std::vector<char>& isColToLeft,
      const MonomialMap<ColIndex>& map
    ) {
      const auto& ring = map.ring();

      // Sort columns by monomial while keeping track of original index.
      MonomialMap<ColIndex>::Reader reader(map);
      typedef std::pair<ColIndex, ConstMonomial> IndexMono;
      std::vector<IndexMono> columns(reader.begin(), reader.end());
      const auto cmp = [&ring](const IndexMono& a, const IndexMono b) {
        return ring.monomialLT(b.second, a.second);
      };
      tbb::parallel_sort(columns.begin(), columns.end(), cmp);

      // Copy monomials and construct projection mapping.
      MATHICGB_ASSERT
        (isColToLeft.size() <= std::numeric_limits<ColIndex>::max());
      ColIndex colCount = static_cast<ColIndex>(isColToLeft.size());
      mProject.resize(isColToLeft.size());
      for (size_t i = 0; i < colCount; ++i) {
        const auto indexMono = columns[i];
        monomial mono = ring.allocMonomial();
        ring.monomialCopy(indexMono.second, mono);

        auto& projected = mProject[indexMono.first];
        projected.left = isColToLeft[indexMono.first];
        if (projected.left) {
          projected.index = static_cast<ColIndex>(mLeftMonomials.size());
          mLeftMonomials.push_back(mono);
        } else {
          projected.index = static_cast<ColIndex>(mRightMonomials.size());
          mRightMonomials.push_back(mono);
        }
      }
      MATHICGB_ASSERT
        (mLeftMonomials.size() + mRightMonomials.size() == isColToLeft.size());
    }

    struct Projected {
      ColIndex index;
      bool left;
    };

    Projected project(const ColIndex index) const {
      MATHICGB_ASSERT(index < mProject.size());
      return mProject[index];
    }

    void project(
      const std::vector<F4ProtoMatrix*>& preBlocks,
      SparseMatrix& left,
      SparseMatrix& right,
      const PolyRing& ring
    ) const {
      left.clear();
      right.clear();
      const auto modulus = static_cast<SparseMatrix::Scalar>(ring.charac());

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
              const auto translated = project(index);
              if (translated.left)
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
              const auto translated = project(index);
              if (translated.left)
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

    void project(
      const std::vector<std::pair<SparseMatrix::Scalar, F4ProtoMatrix::Row>>& from,
      SparseMatrix& left,
      SparseMatrix& right,
      const PolyRing& ring
    ) const {
      left.clear();
      right.clear();
      const auto fromEnd = from.end();
      for (auto fromIt = from.begin(); fromIt != fromEnd; ++fromIt) {
        const auto row = fromIt->second;
        MATHICGB_ASSERT(row.entryCount != 0);
        MATHICGB_ASSERT(row.scalars == 0 || row.externalScalars == 0);
        const auto modulus = static_cast<SparseMatrix::Scalar>(ring.charac());

        if (row.externalScalars != 0) {
          auto indices = row.indices;
          auto indicesEnd = row.indices + row.entryCount;
          auto scalars = row.externalScalars;
          for (; indices != indicesEnd; ++indices, ++scalars) {
            const auto scalar = static_cast<SparseMatrix::Scalar>(*scalars);
            const auto index = *indices;
            const auto translated = project(index);
            if (translated.left)
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
            const auto translated = project(index);
            if (translated.left)
              left.appendEntry(translated.index, *scalars);
            else
              right.appendEntry(translated.index, *scalars);
          }
        }
        const auto rowIndex = left.rowCount();
        MATHICGB_ASSERT(rowIndex == right.rowCount());
        left.rowDone();
        right.rowDone();

        if (fromIt->first != 1) {
          MATHICGB_ASSERT(fromIt->first != 0);
          left.multiplyRow(rowIndex, fromIt->first, modulus);
          right.multiplyRow(rowIndex, fromIt->first, modulus);
          MATHICGB_ASSERT(left.rowBegin(rowIndex).scalar() == 1);
        }

        MATHICGB_ASSERT(left.rowCount() == right.rowCount());
      }
    }

    const std::vector<monomial>& leftMonomials() const {
      return mLeftMonomials;
    }

    std::vector<monomial> moveLeftMonomials() {
      return std::move(mLeftMonomials);
    }

    std::vector<monomial> moveRightMonomials() {
      return std::move(mRightMonomials);
    }

  private:
    std::vector<Projected> mProject;
    std::vector<monomial> mLeftMonomials;
    std::vector<monomial> mRightMonomials;
  };

  class TopBottomProjectionLate {
  public:
    TopBottomProjectionLate(
      const SparseMatrix& left,
      const SparseMatrix& right,
      const SparseMatrix::ColIndex leftColCount,
      const PolyRing& ring
    ):
      mModulus(static_cast<SparseMatrix::Scalar>(ring.charac()))
    {
      const auto noRow = std::numeric_limits<SparseMatrix::ColIndex>::max();
      mTopRows.resize(leftColCount, std::make_pair(0, noRow));

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
          mBottomRows.push_back(std::make_pair(1, row)); //can't be top/reducer
          continue;
        }
        const auto lead = left.rowBegin(row).index();
        if (mTopRows[lead].second != noRow && topEntryCounts[lead]<entryCount)
          mBottomRows.push_back(std::make_pair(1, row)); //other reducer better
        else {
          if (mTopRows[lead].second != noRow)
            mBottomRows.push_back(std::make_pair(1, mTopRows[lead].second));
          topEntryCounts[lead] = entryCount;
          mTopRows[lead].second = row; // do scalar .first later
        }
      }

      const auto modulus = static_cast<SparseMatrix::Scalar>(ring.charac());
      for (SparseMatrix::RowIndex r = 0; r < leftColCount; ++r) {
        const auto row = mTopRows[r].second;
        MATHICGB_ASSERT(row != noRow);
        MATHICGB_ASSERT(left.entryCountInRow(row) > 0);
        MATHICGB_ASSERT(left.rowBegin(row).index() == r);
        MATHICGB_ASSERT(left.rowBegin(row).scalar() != 0);
        MATHICGB_ASSERT(topEntryCounts[r] ==
          left.entryCountInRow(row) + right.entryCountInRow(row));

        const auto leadScalar = left.rowBegin(row).scalar();
        mTopRows[r].first = leadScalar == 1 ? 1 : // 1 is the common case
          modularInverse(leadScalar, modulus);
      }

#ifdef MATHICGB_DEBUG
      for (SparseMatrix::RowIndex r = 0; r < mBottomRows.size(); ++r) {
        const auto row = mBottomRows[r].second;
        MATHICGB_ASSERT(
          left.entryCountInRow(row) + right.entryCountInRow(row) > 0);
        MATHICGB_ASSERT(mBottomRows[r].first == 1);
      }
#endif
    }

    void project(
      SparseMatrix&& in,
      SparseMatrix& top,
      SparseMatrix& bottom
    ) {
      top.clear();
      bottom.clear();

      const auto rowCountTop =
        static_cast<SparseMatrix::RowIndex>(mTopRows.size());
      for (SparseMatrix::RowIndex toRow = 0; toRow < rowCountTop; ++toRow) {
        top.appendRow(in, mTopRows[toRow].second);
        if (mTopRows[toRow].first != 1)
          top.multiplyRow(toRow, mTopRows[toRow].first, mModulus);
      }

      const auto rowCountBottom =
        static_cast<SparseMatrix::RowIndex>(mBottomRows.size());
      for (SparseMatrix::RowIndex toRow = 0; toRow < rowCountBottom; ++toRow) {
        bottom.appendRow(in, mBottomRows[toRow].second);
        if (mBottomRows[toRow].first != 1)
            bottom.multiplyRow(toRow, mBottomRows[toRow].first, mModulus);
      }

      in.clear();
    }

  private:
    const SparseMatrix::Scalar mModulus;
    std::vector<std::pair<SparseMatrix::Scalar, SparseMatrix::RowIndex>> mTopRows;
    std::vector<std::pair<SparseMatrix::Scalar, SparseMatrix::RowIndex>> mBottomRows;
  };

  class TopBottomProjection {
  public:
    TopBottomProjection(
      const std::vector<F4ProtoMatrix*>& blocks,
      const LeftRightProjection& leftRight,
      const PolyRing& ring
    ):
      mReducerRows(leftRight.leftMonomials().size())
    {
      typedef SparseMatrix::RowIndex RowIndex;
      const auto noReducer = std::numeric_limits<RowIndex>::max();
      F4ProtoMatrix::Row noRow = {};
      noRow.indices = 0;
      const auto noCol = std::numeric_limits<SparseMatrix::ColIndex>::max();

      const auto modulus = static_cast<SparseMatrix::Scalar>(ring.charac());

      const auto end = blocks.end();
      for (auto it = blocks.begin(); it != end; ++it) {
        auto& block = **it;
        const auto rowCount = block.rowCount();
        for (RowIndex r = 0; r < rowCount; ++r) {
          const auto row = block.row(r);
          if (row.entryCount == 0)
            continue;

          // Determine leading (minimum index) left entry.
          const auto lead = [&] {
            for (SparseMatrix::ColIndex col = 0; col < row.entryCount; ++col) {
              auto const translated = leftRight.project(row.indices[col]);
              if (translated.left)
                return std::make_pair(col, translated.index);
            }
            return std::make_pair(noCol, noCol);
          }();
          // Decide if this should be a reducer or reducee row.
          if (lead.second == noCol) {
            mReduceeRows.push_back(std::make_pair(1, row)); // no left entries
            continue;
          }
          MATHICGB_ASSERT(row.scalars != 0 || row.externalScalars != 0);

          const auto reducer = mReducerRows[lead.second].second;
          if (
            reducer.entryCount != 0 && // already have a reducer and...
            reducer.entryCount < row.entryCount // ...it is sparser/better
          )
            mReduceeRows.push_back(std::make_pair(1, row));
          else {
            if (reducer.entryCount != 0)
              mReduceeRows.push_back(std::make_pair(1, reducer));
            const auto leadScalar = row.scalars != 0 ? row.scalars[lead.first] :
              static_cast<SparseMatrix::Scalar>
                (row.externalScalars[lead.first]);
            MATHICGB_ASSERT(leadScalar != 0);
            const auto inverse = leadScalar == 1 ?
              1 : modularInverse(leadScalar, modulus);
            mReducerRows[lead.second] = std::make_pair(inverse, row);
          }
        }
      }

#ifdef MATHICGB_DEBUG
  for (size_t  i = 0; i < mReducerRows.size(); ++i) {
    const auto row = mReducerRows[i];
    MATHICGB_ASSERT(row.second.entryCount > 0);
    for (SparseMatrix::ColIndex col = 0; ; ++col) {
      MATHICGB_ASSERT(col < row.second.entryCount);
      const auto projected = leftRight.project(row.second.indices[col]);
      if (projected.left) {
        MATHICGB_ASSERT(projected.index == i);
        const auto leadScalar = row.second.scalars != 0 ?
          row.second.scalars[col] :
          static_cast<SparseMatrix::Scalar>(row.second.externalScalars[col]);
        MATHICGB_ASSERT(modularProduct(leadScalar, row.first, modulus) == 1);
        break;
      }
    }
  }
  for (size_t  i = 0; i < mReduceeRows.size(); ++i) {
    const auto row = mReduceeRows[i];
    MATHICGB_ASSERT(row.second.entryCount > 0);
    MATHICGB_ASSERT(row.first == 1);
  }
#endif
    }

    const std::vector<std::pair<SparseMatrix::Scalar, F4ProtoMatrix::Row>>& reducerRows() const {
      return mReducerRows;
    }

    const std::vector<std::pair<SparseMatrix::Scalar, F4ProtoMatrix::Row>>& reduceeRows() const {
      return mReduceeRows;
    }

  private:
    std::vector<std::pair<SparseMatrix::Scalar, F4ProtoMatrix::Row>> mReducerRows;
    std::vector<std::pair<SparseMatrix::Scalar, F4ProtoMatrix::Row>> mReduceeRows;
  };
}


QuadMatrix F4MatrixProjection::project(
  std::vector<char>&& isColumnToLeft,
  MonomialMap<ColIndex>&& map,
  std::vector<F4ProtoMatrix*>&& matrix
) {
  const auto& ring = map.ring();
  QuadMatrix quadMatrix;
  // Create projections
  LeftRightProjection projection(isColumnToLeft, map);

  if (true) {
    TopBottomProjection topBottom(matrix, projection, ring);

    // Project the pre-blocks into the matrix
    projection.project(topBottom.reducerRows(), quadMatrix.topLeft, quadMatrix.topRight, ring);
    projection.project(topBottom.reduceeRows(), quadMatrix.bottomLeft, quadMatrix.bottomRight, ring);
  } else {
    SparseMatrix left;
    SparseMatrix right;
    projection.project(matrix, left, right, ring);
    TopBottomProjectionLate topBottom(left, right, left.computeColCount(), ring);
    topBottom.project(std::move(left), quadMatrix.topLeft, quadMatrix.bottomLeft);
    topBottom.project(std::move(right), quadMatrix.topRight, quadMatrix.bottomRight);
  }

  quadMatrix.ring = &ring;
  quadMatrix.leftColumnMonomials = projection.moveLeftMonomials();
  quadMatrix.rightColumnMonomials = projection.moveRightMonomials();
  return std::move(quadMatrix);
}
