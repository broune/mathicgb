#include "stdinc.h"
#include "F4MatrixReducer.hpp"

#include "QuadMatrix.hpp"
#include "SparseMatrix.hpp"
#include "PolyRing.hpp"
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <map>
#include <string>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

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

    void addRow(SparseMatrix const& matrix, SparseMatrix::RowIndex row) {
      MATHICGB_ASSERT(row < matrix.rowCount());
      MATHICGB_ASSERT(matrix.colCount() == colCount());
      const auto end = matrix.rowEnd(row);
      for (auto it = matrix.rowBegin(row); it != end; ++it) {
        MATHICGB_ASSERT(it.index() < colCount());
        mEntries[it.index()] = it.scalar();
      }
    }

    template<class Iter>
    void addRowMultiple(
      const SparseMatrix::Scalar multiple,
      const Iter& begin,
      const Iter& end
    ) {
      // Now, you may be wondering why begin and end are passed by reference
      // instead of by value, and that would be a good question. As it turns
      // out, this method does not work otherwise when run in parallel using
      // OpenMP on MS Visual Studio 2012 when being called from reduce().
      // Strange but true.
      //
      // Why is that? To the best of my ability
      // to determine what was going on, it appears that, sometimes,
      // two threads would be running on the same stack when calling this method,
      // overwriting each other's local variables causing all kinds of havoc.
      // My evidence for this is that I could find no other explanation after
      // hours of investigation and that I could consistently get two threads
      // with different return values of omp_get_thread_num() to print out the
      // same address for a local variable - and this would happen just before
      // things went wrong. So at this point I'm concluding that it is a compiler
      // bug. All the writes and reads outside critical sections are to local
      // variables, memory allocated by the same thread or to data structures
      // that do not change within the scope of the parallel code in reduce(),
      // so I don't know what the issue would otherwise be. I thought perhaps
      // not building all the code with OpenMP enabled could be the issue,
      // but changing that did not fix the issue.
      // 
      // Now you may be wondering what that has to do with passing iterators by
      // reference. As it happens, the issue does not happen this way. "But that
      // doesn't make any sense", you say, and you would be right. Feel free
      // to come up with a better explanation of this issue.
      //
      // If you want to take a look at this issue, the issue only turns up for 64
      // bit debug builds. This was on Visual Studio version
      // "11.0.50727.1 RTMREL" - Bjarke Hammersholt Roune


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
      MATHICGB_ASSERT(pivotRow < matrix.rowCount());
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
    SparseMatrix::Scalar modulus,
    const int threadCount
  ) {
    SparseMatrix const& toReduceLeft = qm.bottomLeft;
    SparseMatrix const& toReduceRight = qm.bottomRight;
    SparseMatrix const& reduceByLeft = qm.topLeft;
    SparseMatrix const& reduceByRight = qm.topRight;

    MATHICGB_ASSERT(reduceByLeft.colCount() == reduceByLeft.rowCount());
    const auto pivotCount = reduceByLeft.colCount();
    const auto rowCount = toReduceLeft.rowCount();
    const auto colCountLeft = toReduceLeft.colCount();
    const auto colCountRight = toReduceRight.colCount();

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

    SparseMatrix reduced(colCountRight);

#ifdef _OPENMP
    std::vector<DenseRow<uint64> > denseRowPerThread(threadCount);
#else
    DenseRow<uint64> denseRow;
#endif

    SparseMatrix tmp(pivotCount);

    std::vector<SparseMatrix::RowIndex> rowOrder(rowCount);

#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
    for (OMPIndex rowOMP = 0;
      rowOMP < static_cast<OMPIndex>(rowCount); ++rowOMP) {
      const size_t row = static_cast<size_t>(rowOMP);
#ifdef _OPENMP
      auto& denseRow = denseRowPerThread[omp_get_thread_num()];
#endif
      denseRow.clear(colCountLeft);
      denseRow.addRow(toReduceLeft, row);
      MATHICGB_ASSERT(colCountLeft == pivotCount);

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
#pragma omp critical
      {
        for (size_t pivot = 0; pivot < pivotCount; ++pivot) {
		  MATHICGB_ASSERT(denseRow[pivot] < std::numeric_limits<SparseMatrix::Scalar>::max());
          if (denseRow[pivot] != 0)
            tmp.appendEntry(rowThatReducesCol[pivot], static_cast<SparseMatrix::Scalar>(denseRow[pivot]));
	    }
        tmp.rowDone();
        rowOrder[tmp.rowCount() - 1] = row;
      }
    }

#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
    for (OMPIndex iOMP = 0; iOMP < static_cast<OMPIndex>(rowCount); ++iOMP) {
      const size_t i = static_cast<size_t>(iOMP);
#ifdef _OPENMP
      auto& denseRow = denseRowPerThread[omp_get_thread_num()];
#endif
      size_t row = rowOrder[i];

      denseRow.clear(colCountRight);
      denseRow.addRow(toReduceRight, row);
      auto it = tmp.rowBegin(i);
      const auto end = tmp.rowEnd(i);
      for (; it != end; ++it) {
        const auto begin = reduceByRight.rowBegin(it.index());
        const auto end = reduceByRight.rowEnd(it.index());
        denseRow.addRowMultiple(it.scalar(), begin, end);
      }

#pragma omp critical
      {
        bool zero = true;
	    for (SparseMatrix::ColIndex col = 0; col < colCountRight; ++col) {
          const auto entry =
            static_cast<SparseMatrix::Scalar>(denseRow[col] % modulus);
          if (entry != 0) {
            reduced.appendEntry(col, entry);
            zero = false;
          }
        }
        if (!zero)
          reduced.rowDone();
      }
    }
    return std::move(reduced);
  }

  void reduceToEchelonForm
  (SparseMatrix& toReduce, SparseMatrix::Scalar modulus, int threadCount) {
    // making no assumptions on toReduce except no zero rows

    SparseMatrix::RowIndex const rowCount = toReduce.rowCount();
    SparseMatrix::ColIndex const colCount = toReduce.colCount();

    // dense representation 
    std::vector<DenseRow<uint64> > dense(rowCount);
#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
    for (OMPIndex rowOMP = 0;
      rowOMP < static_cast<OMPIndex>(rowCount); ++rowOMP) {
      const size_t row = static_cast<size_t>(rowOMP);
      MATHICGB_ASSERT(!toReduce.emptyRow(row));
      dense[row].clear(colCount);
      dense[row].addRow(toReduce, row);
    }

    // invariant: all columns in row to the left of leadCols[row] are zero.
    std::vector<SparseMatrix::ColIndex> leadCols(rowCount);

    // pivot rows get copied here before being used to reduce the matrix.
    SparseMatrix reduced;
    reduced.clear(colCount);

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
#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
      for (OMPIndex rowOMP = 0;
        rowOMP < static_cast<OMPIndex>(rowCount); ++rowOMP) {
        const size_t row = static_cast<size_t>(rowOMP);
        MATHICGB_ASSERT(leadCols[row] <= colCount);
        DenseRow<uint64>& denseRow = dense[row];
        if (denseRow.empty())
          continue;

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
#pragma omp critical
          {
            if (!columnHasPivot[col]) {
              columnHasPivot[col] = true;
              isNewReducer = true;
              nextReducers.push_back(std::make_pair(col, row));
            }
          }
          if (isNewReducer)
            denseRow.makeUnitary(modulus, col);
        }
      }
      //std::cout << "done reducing that batch" << std::endl;

      reduced.clear(colCount);
      std::sort(nextReducers.begin(), nextReducers.end());
      for (size_t i = 0; i < nextReducers.size(); ++i) {
        MATHICGB_ASSERT(reduced.colCount() == colCount);
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

#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
    for (OMPIndex rowOMP = 0;
      rowOMP < static_cast<OMPIndex>(rowCount); ++rowOMP) {
      const size_t row = static_cast<size_t>(rowOMP);
      dense[row].takeModulus(modulus);
    }

    toReduce.clear(colCount);
    for (size_t row = 0; row < rowCount; ++row)
      if (!dense[row].empty())
        dense[row].appendTo(toReduce);
  }
}

SparseMatrix F4MatrixReducer::reduce(const QuadMatrix& matrix) {
  MATHICGB_ASSERT(mThreadCount >= 1);
  MATHICGB_ASSERT(matrix.debugAssertValid());
  if (tracingLevel >= 3)
    matrix.printSizes(std::cerr);

  SparseMatrix newPivots(::reduce(matrix, mModulus, mThreadCount));
  ::reduceToEchelonForm(newPivots, mModulus, mThreadCount);
  return std::move(newPivots);
}

namespace {
  /// this has to be a separate function that returns the scalar since signed
  /// overflow is undefine behavior so we cannot check after the cast and
  /// we also cannot set mCharac inside the constructor since it is const.
  SparseMatrix::Scalar checkModulus(const PolyRing& ring) {
    // this assert has to be NO_ASSUME as otherwise the branch below will get
    // optimized out.
    MATHICGB_ASSERT_NO_ASSUME(ring.charac() <=
      std::numeric_limits<SparseMatrix::Scalar>::max());
    if (ring.charac() > std::numeric_limits<SparseMatrix::Scalar>::max())
      throw std::overflow_error("Too large modulus in F4 matrix computation.");
    return static_cast<SparseMatrix::Scalar>(ring.charac());
  }
}

F4MatrixReducer::F4MatrixReducer(const PolyRing& ring, const int threadCount):
  mModulus(checkModulus(ring)),
  mThreadCount(std::max(threadCount, 1)
) {
}
