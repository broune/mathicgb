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

extern int tracingLevel;

template<class T>
class DenseRow {
public:
  DenseRow() {}
  DenseRow(size_t colCount): mEntries(colCount) {}

  /// returns false if all entries are zero
  bool takeModulus(SparseMatrix::Scalar modulus, size_t startCol = 0) {
    typedef typename std::vector<T>::iterator Iter;
    T bitwiseOr = 0;
    Iter end = mEntries.end();
    for (Iter it = mEntries.begin() + startCol; it != end; ++it) {
      if (*it >= modulus)
        *it %= modulus;
      bitwiseOr |= *it;
    }
    return bitwiseOr != 0;
  }

  void clear() {
    mEntries.clear();
  }

  bool empty() const {
    return mEntries.empty();
  }

  void reset(size_t colCount) {
    mEntries.clear();
    mEntries.resize(colCount);
  }

  template<class Iter>
  void addSparseNoModulus(Iter begin, Iter end) {
    for (; begin != end; ++begin) {
      MATHICGB_ASSERT(begin.index() < colCount());
      mEntries[begin.index()] += begin.scalar();
    }
  }

  T& operator[](size_t col) {
    MATHICGB_ASSERT(col < colCount());
    return mEntries[col];
  }
  T const& operator[](size_t col) const {
    MATHICGB_ASSERT(col < colCount());
    return mEntries[col];
  }
  size_t colCount() const {
    return mEntries.size();
  }

  void appendToWithModulus(SparseMatrix& matrix, SparseMatrix::Scalar modulus) {
    matrix.appendRowWithModulus(mEntries, modulus);
  }

  void appendTo(SparseMatrix& matrix, size_t leadCol = 0) {
    matrix.appendRow(mEntries, leadCol);
  }

  void normalize(SparseMatrix::Scalar modulus, size_t lead) {
    MATHICGB_ASSERT(lead < colCount());
    MATHICGB_ASSERT(mEntries[lead] != 0);

    typedef typename std::vector<T>::iterator Iter;
    Iter end = mEntries.end();
    Iter it = mEntries.begin() + lead;
    SparseMatrix::Scalar toInvert = static_cast<SparseMatrix::Scalar>(*it % modulus);
    SparseMatrix::Scalar multiply = modularInverse(toInvert, modulus);
    *it = 1;
    for (++it; it != end; ++it) {
      SparseMatrix::Scalar entry = static_cast<SparseMatrix::Scalar>(*it % modulus);
      if (entry != 0)
        entry = modularProduct(entry, multiply, modulus);
      *it = entry;
    }
  }

  void addRow(SparseMatrix const& matrix, SparseMatrix::RowIndex row) {
    MATHICGB_ASSERT(row < matrix.rowCount());
    MATHICGB_ASSERT(matrix.colCount() == colCount());
    SparseMatrix::RowIterator end = matrix.rowEnd(row);
    for (SparseMatrix::RowIterator it = matrix.rowBegin(row); it != end; ++it) {
      MATHICGB_ASSERT(it.index() < colCount());
      mEntries[it.index()] = it.scalar();
    }
  }

  void addRowPrefix(SparseMatrix const& matrix, SparseMatrix::RowIndex row, size_t stopAtCol) {
    MATHICGB_ASSERT(row < matrix.rowCount());
    MATHICGB_ASSERT(matrix.colCount() == colCount());
    SparseMatrix::RowIterator end = matrix.rowEnd(row);
    for (SparseMatrix::RowIterator it = matrix.rowBegin(row); it != end; ++it) {
      if (it.index() >= stopAtCol)
        break;
      MATHICGB_ASSERT(it.index() < colCount());
      mEntries[it.index()] = it.scalar();
    }
  }

  template<class Iter>
  void addRowMultiple(SparseMatrix::Scalar multiple, Iter begin, Iter end) {
    T mult = static_cast<T>(multiple);
    for (; begin != end; ++begin) {
      MATHICGB_ASSERT(begin.index() < colCount());
      // Watch out for overflow here! This is only OK because
      // T is promoting the multiplication to type T.
      mEntries[begin.index()] += begin.scalar() * mult;
    }
  }

  template<class Iter>
  void addRowPrefixMultiple(SparseMatrix::Scalar multiple, Iter begin, Iter end, SparseMatrix::RowIndex stopAtCol) {
    T mult = static_cast<T>(multiple);
    for (; begin != end; ++begin) {
      if (begin.index() >= stopAtCol)
        break;
      MATHICGB_ASSERT(begin.index() < colCount());
      // Watch out for overflow here! This is only OK because
      // T is promoting the multiplication to type T.
      mEntries[begin.index()] += begin.scalar() * mult;
    }
  }

  void addRowMultiple(SparseMatrix::Scalar multiple, DenseRow const& dense, size_t lead) {
    MATHICGB_ASSERT(dense.colCount() == colCount());
    MATHICGB_ASSERT(lead < colCount());
    T mult = static_cast<T>(multiple);
    size_t colCo = colCount();
    for (size_t col = lead; col < colCo; ++col) {
      // Watch out for overflow here! This is only OK because
      // T is promoting the multiplication to type T.
      mEntries[col] += dense[col] * mult;
    }
  }

  void rowReduceByUnitary
  (size_t pivotRow, SparseMatrix const& matrix, SparseMatrix::Scalar modulus) {
    MATHICGB_ASSERT(pivotRow < matrix.rowCount());
    MATHICGB_ASSERT(matrix.rowBegin(pivotRow).scalar() == 1); // unitary
    MATHICGB_ASSERT(modulus > 1);

    SparseMatrix::RowIterator begin = matrix.rowBegin(pivotRow);
    SparseMatrix::ColIndex col = begin.index();
    SparseMatrix::Scalar entry = mEntries[col] % modulus;
    mEntries[col] = 0;
    if (entry == 0)
      return;
    ++begin; // can skip first entry as we just set it to zero.
    addRowMultiple(modulus - entry, begin, matrix.rowEnd(pivotRow));
  }

  void rowReduceByUnitary
  (size_t pivotRow, const SparseMatrix& matrixLeft, const SparseMatrix& matrixRight, SparseMatrix::Scalar modulus) {
    MATHICGB_ASSERT(pivotRow < matrixLeft.rowCount());
    MATHICGB_ASSERT(pivotRow < matrixRight.rowCount());
    MATHICGB_ASSERT(matrixLeft.rowBegin(pivotRow).scalar() == 1); // unitary
    MATHICGB_ASSERT(modulus > 1);

    SparseMatrix::RowIterator begin = matrixLeft.rowBegin(pivotRow);
    SparseMatrix::ColIndex col = begin.index();
    SparseMatrix::Scalar entry = mEntries[col] % modulus;
    mEntries[col] = 0;
    if (entry == 0)
      return;
    ++begin; // can skip first entry as we just set it to zero.
    addRowMultiple(modulus - entry, begin, matrixLeft.rowEnd(pivotRow)); 

    T mult = modulus - entry;
    auto leftOffset = matrixLeft.colCount();
    MATHICGB_ASSERT(leftOffset <= colCount());
    auto end = matrixRight.rowEnd(pivotRow);
    for (auto it = matrixRight.rowBegin(pivotRow); it != end; ++it) {
      MATHICGB_ASSERT(it.index() < colCount() - leftOffset);
      // Watch out for overflow here! This is only OK because mult is
      // T and so is promoting the multiplication to type T.
      mEntries[it.index() + leftOffset] += it.scalar() * mult;
    }
  }

  void rowReduceByUnitaryPrefix
  (size_t pivotRow, SparseMatrix const& matrix, SparseMatrix::Scalar modulus, SparseMatrix::RowIndex stopAtCol) {
    MATHICGB_ASSERT(pivotRow < matrix.rowCount());
    MATHICGB_ASSERT(matrix.rowBegin(pivotRow).scalar() == 1); // unitary
    MATHICGB_ASSERT(modulus > 1);

    SparseMatrix::RowIterator begin = matrix.rowBegin(pivotRow);
    SparseMatrix::ColIndex col = begin.index();
    if (col >= stopAtCol)
      return;
    SparseMatrix::Scalar entry = mEntries[col] % modulus;
    mEntries[col] = 0;
    if (entry == 0)
      return;
    ++begin; // can skip first entry as we just set it to zero.
    addRowPrefixMultiple(modulus - entry, begin, matrix.rowEnd(pivotRow), stopAtCol); 
  }

  void rowReduceByUnitary(DenseRow const& row, size_t lead, SparseMatrix::Scalar modulus) {
    MATHICGB_ASSERT(row.colCount() == colCount());
    MATHICGB_ASSERT(lead < row.colCount());
    MATHICGB_ASSERT(row[lead] == 1);
    MATHICGB_ASSERT(modulus > 1);

    SparseMatrix::Scalar entry = mEntries[lead] % modulus;
    mEntries[lead] = 0;
    if (entry == 0)
      return;
    addRowMultiple(modulus - entry, row, lead + 1);
  }

  void rowReduce(size_t pivotRow, SparseMatrix const& matrix, SparseMatrix::Scalar modulus) {
    MATHICGB_ASSERT(pivotRow < matrix.rowCount());
    MATHICGB_ASSERT(matrix.rowBegin(pivotRow) != matrix.rowEnd(pivotRow));
    MATHICGB_ASSERT(matrix.rowBegin(pivotRow).scalar() != 0);
    MATHICGB_ASSERT(modulus > 1);

    SparseMatrix::RowIterator begin = matrix.rowBegin(pivotRow);
    SparseMatrix::ColIndex col = begin.index();
    SparseMatrix::Scalar entry = mEntries[col] % modulus;
    mEntries[col] = 0;
    if (entry == 0)
      return;
    SparseMatrix::Scalar reducerLead = begin.scalar();
    ++begin; // can skip first entry as we just set it to zero.
    SparseMatrix::RowIterator end = matrix.rowEnd(pivotRow);
    if (begin == end)
      return;

    SparseMatrix::Scalar inverse = modularInverse(reducerLead, modulus);
    SparseMatrix::Scalar prod = modularProduct(inverse, entry, modulus);
    SparseMatrix::Scalar negProd = modularNegativeNonZero(prod, modulus);
    addRowMultiple(negProd, begin, end);
  }

  // private:
  std::vector<T> mEntries;
};

template<typename Matrix>
void reformMatrix(const Matrix& matA, const Matrix& matB, SparseMatrix& matAB) {
  MATHICGB_ASSERT(matA.rowdim() == matB.rowdim());

  matAB.clear(matA.coldim() + matB.coldim());
  MATHICGB_ASSERT(matAB.colCount() == matA.coldim() + matB.coldim());
  size_t const colCountA = matA.coldim();
  size_t const rowCount = matA.rowdim();

  typedef typename Matrix::Row::const_iterator CIter;
  for (size_t row = 0; row < rowCount; ++row) {
    {
      CIter const endA = matA[row].end();
      for (CIter it = matA[row].begin(); it != endA; ++it) {
        MATHICGB_ASSERT(it->first < colCountA);
        matAB.appendEntry(it->first, it->second);
      }
    }
    {
      CIter const endB = matB[row].end();
      for (CIter it = matB[row].begin(); it != endB; ++it) {
        MATHICGB_ASSERT(it->first < matB.coldim());
        MATHICGB_ASSERT(it->first + colCountA < matAB.colCount());
        matAB.appendEntry(it->first + colCountA, it->second);
      }
    }
    matAB.rowDone();
  }
  MATHICGB_ASSERT(matAB.rowCount() == matA.rowdim());
  MATHICGB_ASSERT(matAB.colCount() == matA.coldim() + matB.coldim());
}

void myReduce
(SparseMatrix const& toReduce,
 SparseMatrix const& reduceBy,
 SparseMatrix::Scalar modulus,
 SparseMatrix& reduced,
 int threadCount) {
  MATHICGB_ASSERT(reduceBy.colCount() >= reduceBy.rowCount());
  MATHICGB_ASSERT(reduceBy.colCount() == toReduce.colCount());
  const auto pivotCount = reduceBy.rowCount();
  const auto colCount = toReduce.colCount();
  const auto rowCount = toReduce.rowCount();

  reduced.clear(toReduce.colCount());

  // pre-calculate what rows are pivots for what columns
  std::vector<SparseMatrix::RowIndex> rowThatReducesCol(pivotCount);
#ifdef MATHICGB_DEBUG
  // fill in an invalid value that can be recognized by asserts to be invalid.
  std::fill(rowThatReducesCol.begin(), rowThatReducesCol.end(), pivotCount);
#endif
  for (SparseMatrix::RowIndex pivot = 0; pivot < pivotCount; ++pivot) {
    MATHICGB_ASSERT(!reduceBy.emptyRow(pivot));
    SparseMatrix::ColIndex col = reduceBy.leadCol(pivot);
    MATHICGB_ASSERT(rowThatReducesCol[col] == pivotCount);
    rowThatReducesCol[col] = pivot;
  }

#ifdef _OPENMP
  std::vector<DenseRow<uint64> > denseRowPerThread(threadCount);
#else
  DenseRow<uint64> denseRow;
#endif

#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (size_t row = 0; row < rowCount; ++row) {
    if (toReduce.emptyRow(row))
      continue;
#ifdef _OPENMP
    DenseRow<uint64>& denseRow = denseRowPerThread[omp_get_thread_num()];
#endif
    denseRow.reset(colCount);
    denseRow.addRow(toReduce, row);
    for (size_t pivot = 0; pivot < pivotCount; ++pivot) {
      if (denseRow[pivot] == 0)
        continue;
      denseRow.rowReduceByUnitary(rowThatReducesCol[pivot], reduceBy, modulus);
    }
    if (denseRow.takeModulus(modulus, pivotCount)) {
#pragma omp critical
      {
        denseRow.appendTo(reduced, pivotCount);
      }
    }
  }
}

void myReduce
(SparseMatrix const& toReduce,
 SparseMatrix const& reduceByLeft,
 SparseMatrix const& reduceByRight,
 SparseMatrix::Scalar modulus,
 SparseMatrix& reduced,
 int threadCount) {
  MATHICGB_ASSERT(reduceByLeft.colCount() == reduceByLeft.rowCount());
  MATHICGB_ASSERT(reduceByLeft.colCount() + reduceByRight.colCount() == toReduce.colCount());
  const auto pivotCount = reduceByLeft.rowCount();
  const auto colCount = toReduce.colCount();
  const auto rowCount = toReduce.rowCount();

  reduced.clear(toReduce.colCount());

  // pre-calculate what rows are pivots for what columns
  std::vector<SparseMatrix::RowIndex> rowThatReducesCol(pivotCount);
#ifdef MATHICGB_DEBUG
  // fill in an invalid value that can be recognized by asserts to be invalid.
  std::fill(rowThatReducesCol.begin(), rowThatReducesCol.end(), pivotCount);
#endif
  for (SparseMatrix::RowIndex pivot = 0; pivot < pivotCount; ++pivot) {
    MATHICGB_ASSERT(!reduceByLeft.emptyRow(pivot));
    SparseMatrix::ColIndex col = reduceByLeft.leadCol(pivot);
    MATHICGB_ASSERT(rowThatReducesCol[col] == pivotCount);
    rowThatReducesCol[col] = pivot;
  }

#ifdef _OPENMP
  std::vector<DenseRow<uint64> > denseRowPerThread(threadCount);
#else
  DenseRow<uint64> denseRow;
#endif

#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (size_t row = 0; row < rowCount; ++row) {
    if (toReduce.emptyRow(row))
      continue;
#ifdef _OPENMP
    DenseRow<uint64>& denseRow = denseRowPerThread[omp_get_thread_num()];
#endif
    denseRow.reset(colCount);
    denseRow.addRow(toReduce, row);
    for (size_t pivot = 0; pivot < pivotCount; ++pivot) {
      if (denseRow[pivot] == 0)
        continue;
      denseRow.rowReduceByUnitary(rowThatReducesCol[pivot], reduceByLeft, reduceByRight, modulus);
    }
    if (denseRow.takeModulus(modulus, pivotCount)) {
#pragma omp critical
      {
        denseRow.appendTo(reduced, pivotCount);
      }
    }
  }
}

void myReduce
(SparseMatrix const& toReduceLeft,
 SparseMatrix const& toReduceRight,
 SparseMatrix const& reduceByLeft,
 SparseMatrix const& reduceByRight,
 SparseMatrix::Scalar modulus,
 SparseMatrix& reduced,
 int threadCount) {
  MATHICGB_ASSERT(reduceByLeft.colCount() == reduceByLeft.rowCount());
  const auto pivotCount = reduceByLeft.rowCount();
  const auto rowCount = toReduceLeft.rowCount();
  const auto colCountLeft = toReduceLeft.colCount();
  const auto colCountRight = toReduceRight.colCount();

  // pre-calculate what rows are pivots for what columns
  std::vector<SparseMatrix::RowIndex> rowThatReducesCol(pivotCount);
#ifdef MATHICGB_DEBUG
  // fill in an invalid value that can be recognized by asserts to be invalid.
  std::fill(rowThatReducesCol.begin(), rowThatReducesCol.end(), pivotCount);
#endif
  for (SparseMatrix::RowIndex pivot = 0; pivot < pivotCount; ++pivot) {
    MATHICGB_ASSERT(!reduceByLeft.emptyRow(pivot));
    SparseMatrix::ColIndex col = reduceByLeft.leadCol(pivot);
    MATHICGB_ASSERT(rowThatReducesCol[col] == pivotCount);
    rowThatReducesCol[col] = pivot;
  }

  reduced.clear(colCountRight);

#ifdef _OPENMP
  std::vector<DenseRow<uint64> > denseRowPerThread(threadCount);
#else
  DenseRow<uint64> denseRow;
#endif

  SparseMatrix tmp;

  std::vector<SparseMatrix::RowIndex> rowOrder(rowCount);

#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (size_t row = 0; row < rowCount; ++row) {
#ifdef _OPENMP
    auto& denseRow = denseRowPerThread[omp_get_thread_num()];
#endif
    denseRow.reset(colCountLeft);
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
      denseRow.addRowMultiple(entry, ++reduceByLeft.rowBegin(row), reduceByLeft.rowEnd(row));
      denseRow[pivot] = entry;
    }
#pragma omp critical
    {
      for (size_t pivot = 0; pivot < pivotCount; ++pivot)
        if (denseRow[pivot] != 0)
          tmp.appendEntry(rowThatReducesCol[pivot], denseRow[pivot]);
      tmp.rowDone();
      rowOrder[tmp.rowCount() - 1] = row;
    }
  }


#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (size_t i = 0; i < rowCount; ++i) {
#ifdef _OPENMP
    auto& denseRow = denseRowPerThread[omp_get_thread_num()];
#endif
    size_t row = rowOrder[i];

    DenseRow<uint64> denseRowX;
    denseRowX.reset(colCountRight);
    denseRowX.addRow(toReduceRight, row);
    auto it = tmp.rowBegin(i);
    auto end = tmp.rowEnd(i);
    for (; it != end; ++it)
      denseRowX.addRowMultiple(it.scalar(), reduceByRight.rowBegin(it.index()), reduceByRight.rowEnd(it.index()));

#pragma omp critical
    {
      bool zero = true;
      for (size_t col = 0; col < colCountRight; ++col) {
        SparseMatrix::Scalar entry = static_cast<SparseMatrix::Scalar>(denseRowX[col] % modulus);
        if (entry != 0)
          reduced.appendEntry(col, entry), zero = false;
      }
      if (!zero)
        reduced.rowDone();
    }
  }


  /*
#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (size_t row = 0; row < rowCount; ++row) {
#ifdef _OPENMP
    auto& denseRow = denseRowPerThread[omp_get_thread_num()];
#endif
    denseRow.reset(colCountRight);
    denseRow.addRow(toReduceRight, row);
    MATHICGB_ASSERT(colCountLeft == pivotCount);
    auto it = tmp.rowBegin(row);
    auto end = tmp.rowEnd(row);
    for (; it != end; ++it)
      denseRow.addRowMultiple(it.scalar(), reduceByRight.rowBegin(it.index()), reduceByRight.rowEnd(it.index()));
#pragma omp critical
    {
      bool zero = true;
      for (size_t col = 0; col < colCountRight; ++col) {
        SparseMatrix::Scalar entry = static_cast<SparseMatrix::Scalar>(denseRow[col] % modulus);
        if (entry != 0)
          reduced.appendEntry(col, entry), zero = false;
      }
      if (!zero)
        reduced.rowDone();
    }
  }
*/
}

void myReduceToEchelonForm5
(SparseMatrix& toReduce, SparseMatrix::Scalar modulus, size_t threadCount) {
  // making no assumptions on toReduce except no zero rows

  SparseMatrix::RowIndex const rowCount = toReduce.rowCount();
  SparseMatrix::ColIndex const colCount = toReduce.colCount();

  // dense representation 
  std::vector<DenseRow<uint64> > dense(rowCount);
#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (SparseMatrix::RowIndex row = 0; row < rowCount; ++row) {
    MATHICGB_ASSERT(!toReduce.emptyRow(row));
    dense[row].reset(colCount);
    dense[row].addRow(toReduce, row);
  }

  // invariant: all columns in row to the left of leadCols[row] are zero.
  std::vector<size_t> leadCols(rowCount);

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
    for (size_t row = 0; row < rowCount; ++row) {
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
      size_t col;
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
          denseRow.normalize(modulus, col);
      }
    }
    //std::cout << "done reducing that batch" << std::endl;

    reduced.clear(colCount);
    std::sort(nextReducers.begin(), nextReducers.end());
    for (size_t i = 0; i < nextReducers.size(); ++i) {
      MATHICGB_ASSERT(reduced.colCount() == colCount);
      size_t const row = nextReducers[i].second;

#ifdef MATHICGB_DEBUG
      size_t const col = nextReducers[i].first;
#endif
      MATHICGB_ASSERT(col == nextReducers[i].first);
      MATHICGB_ASSERT(columnHasPivot[col]);
      MATHICGB_ASSERT(dense[row].colCount() == colCount);
      MATHICGB_ASSERT(dense[row][col] == 1); // already normalized
      MATHICGB_ASSERT(reduced.rowCount() == i);
      MATHICGB_ASSERT(!isPivotRow[row]);

      dense[row].appendTo(reduced); // already nornamlized
      isPivotRow[row] = true;
    }
    nextReducers.clear();
  }

#pragma omp parallel for num_threads(threadCount) schedule(dynamic)
  for (size_t row = 0; row < rowCount; ++row)
    dense[row].takeModulus(modulus);

  toReduce.clear(colCount);
  for (size_t row = 0; row < rowCount; ++row)
    if (!dense[row].empty())
      dense[row].appendTo(toReduce);
}

class IOException : public std::runtime_error {
 public:
 IOException(): std::runtime_error("File error.") {}
};

template<class T>
T readOne(FILE* file) {
  T t;
  if(fread(&t, sizeof(T), 1, file) != 1)
    throw IOException();
  return t;
}

template<class T>
void writeOne(const T& t, FILE* file) {
  if (fwrite(&t, sizeof(T), 1, file) != 1)
    throw IOException();
}

template<class T>
void writeMany(const std::vector<T>& v, FILE* file) {
  if (v.empty())
    return;
  if(fwrite(&v[0], sizeof(T), v.size(), file) != v.size())
    throw IOException();  
}

template<class T>
void readMany(FILE* file, size_t count, std::vector<T>& v) {
  size_t const originalSize = v.size();
  v.resize(originalSize + count);
  if(fread(&v[originalSize], sizeof(T), count, file) != count)
    throw IOException();
}

// Writes an SparseMatrix
void writeSparseMatrix
(const SparseMatrix& matrix, SparseMatrix::Scalar modulus, const std::string& fileName) {
  const uint32 rowCount = matrix.rowCount();
  const uint32 colCount = matrix.colCount();
  const uint64 entryCount = matrix.entryCount();

  FILE* file = fopen(fileName.c_str(), "w");
  if (file == NULL)
    throw IOException();

  writeOne<uint32>(rowCount, file);
  writeOne<uint32>(colCount, file);
  writeOne<uint32>(modulus, file);
  writeOne<uint64>(entryCount, file);

  writeMany<uint16>(matrix.entries(), file);
  writeMany<uint32>(matrix.colIndices(), file);

  std::vector<uint32> sizes;
  for (size_t row = 0; row < rowCount; ++row)
    sizes.push_back(matrix.entryCountInRow(row));
  writeMany<uint32>(sizes, file);

  // todo: don't leak file on exception.
  fclose(file);
}

// Reads an SparseMatrix without any fseeks and returns the modulus.
SparseMatrix::Scalar readSparseMatrix(const std::string& fileName, SparseMatrix& matrix)
{
  FILE* file = fopen(fileName.c_str(), "r");
  if (file == NULL)
    throw IOException();

  uint32 const rowCount = readOne<uint32>(file);
  uint32 const colCount = readOne<uint32>(file);
  uint32 const modulus = readOne<uint32>(file);
  uint64 const entryCount = readOne<uint64>(file);

  std::vector<uint16> entries;
  readMany(file, entryCount, entries);

  std::vector<uint32> indices;
  readMany(file, entryCount, indices);

  std::vector<uint32> sizes;
  readMany(file, rowCount, sizes);

  matrix.setToAndTakeMemory(indices, entries, sizes, colCount);
  return modulus;
}

// doesn't need to be fast.
int integerLog10(unsigned int val) {
  int ret = -1;
  while (val != 0) {
    val /= 10;
    ret++;
  }
  return ret;
}

size_t integerPow10(int l) {
  if (l < 0)
    return 0;
  size_t ret = 1;
  for (; l > 0; --l)
    ret *= 10;
  return ret;
}

void printMap(std::map<int, size_t> const& m, std::ostream& out) {
  std::map<int, size_t>::const_iterator it = m.begin();
  std::map<int, size_t>::const_iterator end = m.end();
  for (; it != end; ++it)
    out << integerPow10(it->first) << ":\t" << it->second << '\n';  
}

void makeLogHistogram
(std::vector<size_t> numbers,
 std::map<int, size_t>& logHis) {
  logHis.clear();
  std::vector<size_t>::const_iterator it = numbers.begin();
  std::vector<size_t>::const_iterator end = numbers.end();
  for (; it != end; ++it)
    ++logHis[integerLog10(*it)];
}

void computeDensities
(SparseMatrix const& matrix,
 std::vector<size_t>& rowDensity,
 std::vector<size_t>& colDensity) {
  rowDensity.clear();
  rowDensity.resize(matrix.rowCount());
  colDensity.clear();
  colDensity.resize(matrix.colCount());

  size_t const rowCount = matrix.rowCount();
  for (size_t row = 0; row < rowCount; ++row) {
    MATHICGB_ASSERT(row < rowDensity.size());
    SparseMatrix::RowIterator it = matrix.rowBegin(row);
    SparseMatrix::RowIterator end = matrix.rowEnd(row);
    for (; it != end; ++it) {
      MATHICGB_ASSERT(it.index() < colDensity.size());
      MATHICGB_ASSERT(it.scalar() != 0);
      ++rowDensity[row];
      ++colDensity[it.index()];
    }
  }
}


void printLogDensityHistograms
(SparseMatrix const& matrix, std::ostream& out, const char* name = 0) {
  std::vector<size_t> rowDensities;
  std::vector<size_t> colDensities;
  computeDensities(matrix, rowDensities, colDensities);

  if (name != 0)
    out << "/********** Density histogram (" << name << ") **********\n";
  else
    out << "/********** Density histogram ***********\n";
      
  out << "Matrix has " << matrix.rowCount()
      << " rows. log2(rows) = " << integerLog10(matrix.rowCount()) << '\n'
      << "Matrix has " << matrix.colCount()
      << " cols. log2(cols) = " << integerLog10(matrix.colCount()) << '\n';

  out << "\nLog row density histogram:\ndens\trows\n";
  std::map<int, size_t> rowLogHis;
  makeLogHistogram(rowDensities, rowLogHis);
  printMap(rowLogHis, out);

  out << "\nLog col density histogram:\ndens\tcol\n";
  std::map<int, size_t> colLogHis;
  makeLogHistogram(colDensities, colLogHis);
  printMap(colLogHis, out);

  out << "\\****************************************\n";
}

// Takes the pivot rows from matrix and copies them into pivots. The
// remaining rows go into nonPivots. Also reorders columns so that the
// left part of pivots is upper triangular.
// Also sorts rows of nonPivots to be densest first.
// Ignores zero rows.
void spliceMatrix(const SparseMatrix& matrix, SparseMatrix& pivots, SparseMatrix& nonPivots) {
  const SparseMatrix::RowIndex rowCount = matrix.rowCount();
  const SparseMatrix::ColIndex colCount = matrix.colCount();

  static const SparseMatrix::RowIndex noPivot =
    std::numeric_limits<SparseMatrix::RowIndex>::max();
  std::vector<SparseMatrix::RowIndex> pivotRowOfCol(colCount);
  std::fill(pivotRowOfCol.begin(), pivotRowOfCol.end(), noPivot);

  // determine pivot rows and columns
  for (size_t row = 0; row < rowCount; ++row) {
    const SparseMatrix::ColIndex entryCount = matrix.entryCountInRow(row);
    if (entryCount == 0)
      continue; // ignore zero rows
    const SparseMatrix::ColIndex pivotCol = matrix.leadCol(row);
    SparseMatrix::RowIndex& pivotRow = pivotRowOfCol[pivotCol];
    if (pivotRow == noPivot || entryCount < matrix.entryCountInRow(pivotRow))
      pivotRow = row; // prefer sparse pivots
  }

  // permutation of columns to put pivots left without reordering
  // columns in any other way.
  std::vector<SparseMatrix::ColIndex> colPerm(colCount);
  SparseMatrix::RowIndex columnsDecided = 0;

  // choice of rows to make left of pivots matrix upper triangular
  std::vector<SparseMatrix::RowIndex> pivotRows;

  // Go through pivot columns to compute perm and pivotRows
  for (size_t col = 0; col < colCount; ++col) {
    if (pivotRowOfCol[col] != noPivot) {
      colPerm[col] = columnsDecided; // pivot columns first
      ++columnsDecided;
      pivotRows.push_back(pivotRowOfCol[col]);
    }
  }
  SparseMatrix::RowIndex minNonPivotCol = columnsDecided;

  for (size_t col = 0; col < colCount; ++col) {
    if (pivotRowOfCol[col] == noPivot) {
      colPerm[col] = columnsDecided; // non-pivot columns last
      ++columnsDecided;
    }
  }
  MATHICGB_ASSERT(columnsDecided == colCount);

  // choice of rows to make pivots matrix sorted by decreasing density.
  std::vector<SparseMatrix::RowIndex> nonPivotRows;

  // put density first in a pair to make sorting easy.
  // TODO: use sorting object instead of making pairs like this.
  std::vector<std::pair<SparseMatrix::ColIndex, SparseMatrix::RowIndex> > nonPivotData;
  for (size_t row = 0; row < rowCount; ++row) {
    const SparseMatrix::ColIndex entryCount = matrix.entryCountInRow(row);
    if (entryCount == 0)
      continue; // ignore zero rows
    if (row != pivotRowOfCol[matrix.leadCol(row)])
      nonPivotData.push_back(std::make_pair(entryCount, row));
  }
  std::sort(nonPivotRows.rbegin(), nonPivotRows.rend());
  for (size_t i = 0; i < nonPivotData.size(); ++i)
    nonPivotRows.push_back(nonPivotData[i].second);

  // create matrices
  struct LocalFunction {
    void makeMatrix(const SparseMatrix& matrix,
                    const std::vector<SparseMatrix::RowIndex>& rowIndices,
                    const std::vector<SparseMatrix::ColIndex>& colPerm,
                    const SparseMatrix::ColIndex minNonPivotCol,
                    SparseMatrix& out) {
      typedef std::vector<SparseMatrix::RowIndex>::const_iterator Iter;
      Iter end = rowIndices.end();

      // reserve space
      out.clear(matrix.colCount());
      size_t entryCount = 0;
      for (Iter it = rowIndices.begin(); it != end; ++it)
        entryCount += matrix.entryCountInRow(*it);
      out.reserveEntries(entryCount);
      out.reserveRows(rowIndices.size());

      for (Iter it = rowIndices.begin(); it != end; ++it) {
        // Do two passes to avoid having to sort indices. They will
        // be in increasing order in this way.
        SparseMatrix::RowIterator begin = matrix.rowBegin(*it);
        SparseMatrix::RowIterator end = matrix.rowEnd(*it);
        for (SparseMatrix::RowIterator it = begin; it != end; ++it)
          if (colPerm[it.index()] < minNonPivotCol) // pivot columns first
            out.appendEntry(colPerm[it.index()], it.scalar());
        for (SparseMatrix::RowIterator it = begin; it != end; ++it)
          if (colPerm[it.index()] >= minNonPivotCol) // then non-pivot columns
            out.appendEntry(colPerm[it.index()], it.scalar());
        out.rowDone();
      }      
    }
  } f; // static function on local structs not allowed :-(
  f.makeMatrix(matrix, pivotRows, colPerm, minNonPivotCol, pivots);
  f.makeMatrix(matrix, nonPivotRows, colPerm, minNonPivotCol, nonPivots);
}

void concatenateMatricesHorizontal
(const SparseMatrix& a, const SparseMatrix& b, SparseMatrix& concatenation) {
  MATHICGB_ASSERT(a.rowCount() == b.rowCount());
  // todo: check overflow of colcount type
  const SparseMatrix::ColIndex bOffset = a.colCount();
  concatenation.clear(a.colCount() + b.colCount());
  if (concatenation.colCount() < a.colCount()) {
    MATHICGB_ASSERT(false);
    throw std::overflow_error
      ("Too many columns in matrices being concatenated.");
  }
  
  const SparseMatrix::ColIndex colBOffset = a.colCount();
  const SparseMatrix::RowIndex rowCount = a.rowCount();
  for (SparseMatrix::RowIndex row = 0; row < rowCount; ++row) {
    {
      auto end = a.rowEnd(row);
      for (auto it = a.rowBegin(row); it != end; ++it)
        concatenation.appendEntry(it.index(), it.scalar());
    }
    {
      auto end = b.rowEnd(row);
      for (auto it = b.rowBegin(row); it != end; ++it)
        concatenation.appendEntry(it.index() + colBOffset, it.scalar());
    }
    concatenation.rowDone();
  }
}

void F4MatrixReducer::reduce
(const PolyRing& ring, QuadMatrix& matrix, SparseMatrix& newPivots) {
  if (tracingLevel >= 3)
    std::cerr << "Row reducing (" << matrix.topLeft.rowCount()
              << " + " << matrix.bottomLeft.rowCount()
              << ") x (" << matrix.topLeft.colCount()
              << " + " << matrix.topRight.colCount()
              << ") matrix.\n#nz entry: top-left "
              << matrix.topLeft.entryCount()
              << ", top right " << matrix.topRight.entryCount()
              << ", bot left " << matrix.bottomLeft.entryCount()
              << ", bot right " << matrix.bottomRight.entryCount()
              << '\n';

  if (ring.charac() > std::numeric_limits<SparseMatrix::Scalar>::max())
    throw std::overflow_error("Too large modulus in F4 matrix computation.");

  SparseMatrix::Scalar modulus = ring.charac();

  {

    const SparseMatrix::ColIndex pivotColCount = matrix.topLeft.colCount();

    //SparseMatrix matrixCD;    concatenateMatricesHorizontal    (matrix.bottomLeft, matrix.bottomRight, matrixCD);

    myReduce(matrix.bottomLeft, matrix.bottomRight, matrix.topLeft, matrix.topRight, modulus, newPivots, mThreadCount);
    //    myReduce(matrixCD, matrix.topLeft, matrix.topRight, modulus, newPivots, mThreadCount);    newPivots.trimLeadingZeroColumns(pivotColCount);
  }

  myReduceToEchelonForm5(newPivots, modulus, mThreadCount);
}
