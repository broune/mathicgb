#include "stdinc.h"
#include "F4MatrixReducer.hpp"

#include "QuadMatrix.hpp"
#include "SparseMatrix.hpp"
#include "PolyRing.hpp"
#include <tbb/tbb.h>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <map>
#include <string>
#include <cstdio>
#include <iostream>

namespace {
  template<class T>
  class DenseRow {
  public:
    DenseRow() {}
    DenseRow(size_t colCount): mEntries(colCount) {}

    /// returns false if all entries are zero
    bool takeModulus(const SparseMatrix::Scalar modulus) {
      T bitwiseOr = 0; // bitwise or of all entries after modulus
      const auto end = mEntries.end();
      for (auto it = mEntries.begin(); it != end; ++it) {
        if (*it >= modulus)
          *it %= modulus;
        bitwiseOr |= *it;
      }
      return bitwiseOr != 0;
    }

    size_t colCount() const {return mEntries.size();}
    bool empty() const {return mEntries.empty();}

    void clear(size_t colCount = 0) {
      mEntries.clear();
      mEntries.resize(colCount);
    }

    T& operator[](size_t col) {
      MATHICGB_ASSERT(col < colCount());
      return mEntries[col];
    }

    T const& operator[](size_t col) const {
      MATHICGB_ASSERT(col < colCount());
      return mEntries[col];
    }

    void appendTo(SparseMatrix& matrix) {matrix.appendRow(mEntries);}

    void makeUnitary(const SparseMatrix::Scalar modulus, const size_t lead) {
      MATHICGB_ASSERT(lead < colCount());
      MATHICGB_ASSERT(mEntries[lead] != 0);

      const auto end = mEntries.end();
      auto it = mEntries.begin() + lead;
      const auto toInvert = static_cast<SparseMatrix::Scalar>(*it % modulus);
      const auto multiply = modularInverse(toInvert, modulus);
      *it = 1;
      for (++it; it != end; ++it) {
        const auto entry = static_cast<SparseMatrix::Scalar>(*it % modulus);
        if (entry != 0)
          *it = modularProduct(entry, multiply, modulus);
        else
          *it = entry;
      }
    }

    void addRow(const SparseMatrix& matrix, SparseMatrix::RowIndex row) {
      MATHICGB_ASSERT(row < matrix.rowCount());
      const auto end = matrix.rowEnd(row);
      for (auto it = matrix.rowBegin(row); it != end; ++it) {
        MATHICGB_ASSERT(it.index() < colCount());
        mEntries[it.index()] = it.scalar();
      }
    }

    template<class Iter>
    void addRowMultiple(
      const SparseMatrix::Scalar multiple,
      const Iter begin,
      const Iter end
    ) {
      // MATHICGB_RESTRICT on entries is important. It fixed a performance
      // regression on MSVC 2012 which otherwise was not able to determine that
      // entries is not an alias for anything else in this loop. I suspect that
      // this is because MSVC 2012 does not make use of the strict aliasing
      // rule. The regression occurred when reusing the DenseRow object instead
      // of making a new one. I suspect that MSVC 2012 was then previously able
      // to tell that entries is not an alias since new does not return
      // aliases.
      T* const MATHICGB_RESTRICT entries = mEntries.data();
      for (Iter it = begin; it != end; ++it) {
        MATHICGB_ASSERT(it.index() < colCount());
        MATHICGB_ASSERT(entries + it.index() == &mEntries[it.index()]);
        entries[it.index()] += it.scalar() * static_cast<T>(multiple);
      }
    }

    void rowReduceByUnitary(
      const size_t pivotRow,
      const SparseMatrix& matrix,
      const SparseMatrix::Scalar modulus
    ) {
      MATHICGB_ASSERT(matrix.rowBegin(pivotRow).scalar() == 1); // unitary
      MATHICGB_ASSERT(modulus > 1);

      auto begin = matrix.rowBegin(pivotRow);
      const SparseMatrix::ColIndex col = begin.index();
      const SparseMatrix::Scalar entry = mEntries[col] % modulus;
      mEntries[col] = 0;
      if (entry == 0)
        return;
      ++begin; // can skip first entry as we just set it to zero.
      addRowMultiple(modulus - entry, begin, matrix.rowEnd(pivotRow));
    }

  private:
    std::vector<T> mEntries;
  };

  SparseMatrix reduce(
    const QuadMatrix& qm,
    SparseMatrix::Scalar modulus
  ) {
    const SparseMatrix& toReduceLeft = qm.bottomLeft;
    const SparseMatrix& toReduceRight = qm.bottomRight;
    const SparseMatrix& reduceByLeft = qm.topLeft;
    const SparseMatrix& reduceByRight = qm.topRight;

    const auto leftColCount = qm.computeLeftColCount();
//      static_cast<SparseMatrix::ColIndex>(qm.leftColumnMonomials.size());
    const auto rightColCount = static_cast<SparseMatrix::ColIndex>(qm.computeRightColCount());
//      static_cast<SparseMatrix::ColIndex>(qm.rightColumnMonomials.size());
    MATHICGB_ASSERT(leftColCount == reduceByLeft.rowCount());
    const auto pivotCount = leftColCount;
    const auto rowCount = toReduceLeft.rowCount();

    // ** pre-calculate what rows are pivots for what columns.

    // Store column indexes instead of row indices as the matrix is square
    // anyway (so all indices fit) and we are going to store this as a column
    // index later on.
    std::vector<SparseMatrix::ColIndex> rowThatReducesCol(pivotCount);
#ifdef MATHICGB_DEBUG
    // fill in an invalid value that can be recognized by asserts to be invalid.
    std::fill(rowThatReducesCol.begin(), rowThatReducesCol.end(), pivotCount);
#endif
    for (SparseMatrix::ColIndex pivot = 0; pivot < pivotCount; ++pivot) {
      MATHICGB_ASSERT(!reduceByLeft.emptyRow(pivot));
      SparseMatrix::ColIndex col = reduceByLeft.leadCol(pivot);
      MATHICGB_ASSERT(rowThatReducesCol[col] == pivotCount);
      rowThatReducesCol[col] = pivot;
    }

    SparseMatrix reduced(qm.topRight.memoryQuantum());

    tbb::enumerable_thread_specific<DenseRow<uint64>> denseRowPerThread([&](){
      return DenseRow<uint64>();
    }); 

    SparseMatrix tmp(qm.topRight.memoryQuantum());

    std::vector<SparseMatrix::RowIndex> rowOrder(rowCount);

    tbb::mutex lock;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, rowCount),
      [&](const tbb::blocked_range<size_t>& range)
      {for (auto it = range.begin(); it != range.end(); ++it)
    {
      const size_t row = it;
      auto& denseRow = denseRowPerThread.local();

      denseRow.clear(leftColCount);
      denseRow.addRow(toReduceLeft, row);
      MATHICGB_ASSERT(leftColCount == pivotCount);

      for (size_t pivot = 0; pivot < pivotCount; ++pivot) {
        if (denseRow[pivot] == 0)
          continue;
        auto entry = denseRow[pivot];
        entry %= modulus;
        if (entry == 0) {
          denseRow[pivot] = 0;
          continue;
        }
        entry = modulus - entry;
        const auto row = rowThatReducesCol[pivot];
        MATHICGB_ASSERT(row < pivotCount);
        MATHICGB_ASSERT(!reduceByLeft.emptyRow(row));
        MATHICGB_ASSERT(reduceByLeft.leadCol(row) == pivot);
        MATHICGB_ASSERT(entry < std::numeric_limits<SparseMatrix::Scalar>::max());
        denseRow.addRowMultiple(static_cast<SparseMatrix::Scalar>(entry), ++reduceByLeft.rowBegin(row), reduceByLeft.rowEnd(row));
        denseRow[pivot] = entry;
      }
      tbb::mutex::scoped_lock lockGuard(lock);
      for (size_t pivot = 0; pivot < pivotCount; ++pivot) {
		MATHICGB_ASSERT(denseRow[pivot] < std::numeric_limits<SparseMatrix::Scalar>::max());
        if (denseRow[pivot] != 0)
          tmp.appendEntry(rowThatReducesCol[pivot], static_cast<SparseMatrix::Scalar>(denseRow[pivot]));
	  }
      tmp.rowDone();
      rowOrder[tmp.rowCount() - 1] = row;
    }});

    tbb::parallel_for(tbb::blocked_range<size_t>(0, rowCount),
      [&](const tbb::blocked_range<size_t>& range)
      {for (auto iter = range.begin(); iter != range.end(); ++iter)
    {
      const size_t i = iter;
      const size_t row = rowOrder[i];
      auto& denseRow = denseRowPerThread.local();

      denseRow.clear(rightColCount);
      denseRow.addRow(toReduceRight, row);
      auto it = tmp.rowBegin(i);
      const auto end = tmp.rowEnd(i);
      for (; it != end; ++it) {
        const auto begin = reduceByRight.rowBegin(it.index());
        const auto end = reduceByRight.rowEnd(it.index());
        denseRow.addRowMultiple(it.scalar(), begin, end);
      }

      tbb::mutex::scoped_lock lockGuard(lock);
      bool zero = true;
	  for (SparseMatrix::ColIndex col = 0; col < rightColCount; ++col) {
        const auto entry =
          static_cast<SparseMatrix::Scalar>(denseRow[col] % modulus);
        if (entry != 0) {
          reduced.appendEntry(col, entry);
          zero = false;
        }
      }
      if (!zero)
        reduced.rowDone();
    }});
    return std::move(reduced);
  }

  SparseMatrix reduceToEchelonForm(
    const SparseMatrix& toReduce,
    const SparseMatrix::Scalar modulus
  ) {
    const auto colCount = toReduce.computeColCount();
    const auto rowCount = toReduce.rowCount();

    // convert to dense representation 
    std::vector<DenseRow<uint64>> dense(rowCount);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, rowCount),
      [&](const tbb::blocked_range<size_t>& range)
      {for (auto it = range.begin(); it != range.end(); ++it)
    {
      const size_t row = it;
      if (toReduce.emptyRow(row))
        return;
      dense[row].clear(colCount);
      dense[row].addRow(toReduce, row);
    }});

    // invariant: all columns in row to the left of leadCols[row] are zero.
    std::vector<SparseMatrix::ColIndex> leadCols(rowCount);

    // pivot rows get copied here before being used to reduce the matrix.
    SparseMatrix reduced(toReduce.memoryQuantum());

    // (col,row) in nextReducers, then use row as a pivot in column col
    // for the next iteration.
    std::vector<std::pair<SparseMatrix::ColIndex, SparseMatrix::RowIndex> > nextReducers;

    // isPivotRow[row] is true if row is or has been used as a pivot.
    std::vector<bool> isPivotRow(rowCount);

    // columnHasPivot[col] is true if a pivot row for column col has
    // been chosen.
    std::vector<bool> columnHasPivot(colCount);

    bool firstIteration = true;
    while (firstIteration || reduced.rowCount() > 0) {
      firstIteration = false;
      size_t const reducerCount = reduced.rowCount();

      //std::cout << "reducing " << reduced.rowCount() << " out of " << toReduce.rowCount() << std::endl;
      tbb::mutex lock;
      tbb::parallel_for(tbb::blocked_range<size_t>(0, rowCount),
        [&](const tbb::blocked_range<size_t>& range)
        {for (auto it = range.begin(); it != range.end(); ++it)
      {
        const size_t row = it;
        MATHICGB_ASSERT(leadCols[row] <= colCount);
        DenseRow<uint64>& denseRow = dense[row];
        if (denseRow.empty())
          return;

        // reduce by each row of reduced.
        for (size_t reducerRow = 0; reducerRow < reducerCount; ++reducerRow) {
          size_t const col = reduced.rowBegin(reducerRow).index();
          if (denseRow[col] == 0 || (isPivotRow[row] && col == leadCols[row]))
            continue;
          denseRow.rowReduceByUnitary(reducerRow, reduced, modulus);
        }

        // update leadCols[row]
        SparseMatrix::ColIndex col;
        MATHICGB_ASSERT(leadCols[row] <= colCount);
        for (col = leadCols[row]; col < colCount; ++col) {
          denseRow[col] %= modulus;
          if (denseRow[col] != 0)
            break;
        }
        leadCols[row] = col;
        MATHICGB_ASSERT(leadCols[row] <= colCount);

        // note if we have found a new pivot row
        if (col == colCount)
          denseRow.clear();
        else {
          MATHICGB_ASSERT(col < colCount);
          bool isNewReducer = false;
          {
            tbb::mutex::scoped_lock lockGuard(lock);
            if (!columnHasPivot[col]) {
              columnHasPivot[col] = true;
              isNewReducer = true;
              nextReducers.push_back(std::make_pair(col, row));
            }
          }
          if (isNewReducer)
            denseRow.makeUnitary(modulus, col);
        }
      }});
      //std::cout << "done reducing that batch" << std::endl;

      reduced.clear();
      std::sort(nextReducers.begin(), nextReducers.end());
      for (size_t i = 0; i < nextReducers.size(); ++i) {
        size_t const row = nextReducers[i].second;

        MATHICGB_ASSERT(static_cast<bool>
          (columnHasPivot[nextReducers[i].first]));
        MATHICGB_ASSERT(dense[row].colCount() == colCount);
        MATHICGB_ASSERT(dense[row][nextReducers[i].first] == 1);
        MATHICGB_ASSERT(reduced.rowCount() == i);
        MATHICGB_ASSERT(!isPivotRow[row]);

        dense[row].appendTo(reduced); // already unitary
        isPivotRow[row] = true;
      }
      nextReducers.clear();
    }

    tbb::parallel_for(tbb::blocked_range<size_t>(0, rowCount),
      [&](const tbb::blocked_range<size_t>& range)
      {for (auto it = range.begin(); it != range.end(); ++it)
    {
      const size_t row = it;
      dense[row].takeModulus(modulus);
    }});

    reduced.clear();
    for (size_t row = 0; row < rowCount; ++row)
      if (!dense[row].empty())
        dense[row].appendTo(reduced);
    return std::move(reduced);
  }
}


SparseMatrix F4MatrixReducer::reduceToBottomRight(const QuadMatrix& matrix) {
  MATHICGB_ASSERT(matrix.debugAssertValid());
  if (tracingLevel >= 3)
    matrix.printSizes(std::cerr);
  return reduce(matrix, mModulus);
}

SparseMatrix F4MatrixReducer::reducedRowEchelonForm(
  const SparseMatrix& matrix
) {
  return reduceToEchelonForm(matrix, mModulus);
}

SparseMatrix F4MatrixReducer::reducedRowEchelonFormBottomRight(
  const QuadMatrix& matrix
) {
  return reducedRowEchelonForm(reduceToBottomRight(matrix));
}

namespace {
  /// this has to be a separate function that returns the scalar since signed
  /// overflow is undefine behavior so we cannot check after the cast and
  /// we also cannot set the modulus field inside the constructor since it is
  /// const.
  SparseMatrix::Scalar checkModulus(const coefficient modulus) {
    // this assert has to be NO_ASSUME as otherwise the branch below will get
    // optimized out.
    MATHICGB_ASSERT_NO_ASSUME(modulus <=
      std::numeric_limits<SparseMatrix::Scalar>::max());
    if (modulus > std::numeric_limits<SparseMatrix::Scalar>::max())
      throw std::overflow_error("Too large modulus in F4 matrix reduction.");
    return static_cast<SparseMatrix::Scalar>(modulus);
  }
}

F4MatrixReducer::F4MatrixReducer(const coefficient modulus):
  mModulus(checkModulus(modulus)) {}
