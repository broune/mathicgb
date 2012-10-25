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

  void appendTo(SparseMatrix& matrix, SparseMatrix::ColIndex leadCol = 0) {
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
    const auto end = matrix.rowEnd(row);
    for (auto it = matrix.rowBegin(row); it != end; ++it) {
      MATHICGB_ASSERT(it.index() < colCount());
      mEntries[it.index()] = it.scalar();
    }
  }

  void addRowPrefix(SparseMatrix const& matrix, SparseMatrix::RowIndex row, size_t stopAtCol) {
    MATHICGB_ASSERT(row < matrix.rowCount());
    MATHICGB_ASSERT(matrix.colCount() == colCount());
    const auto end = matrix.rowEnd(row);
    for (auto it = matrix.rowBegin(row); it != end; ++it) {
      if (it.index() >= stopAtCol)
        break;
      MATHICGB_ASSERT(it.index() < colCount());
      mEntries[it.index()] = it.scalar();
    }
  }


  template<class Iter>
  void addRowMultiple(SparseMatrix::Scalar multiple, const Iter& begin, const Iter& end) {
    // Now, you may be wondering why begin and end are passed by reference
    // instead of by value, and that would be a good question. As it turns
    // out, this method does not work otherwise when run in parallel using
    // OpenMP on MS Visual Studio 2012 when being called from myReduce().
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
    // that do not change within the scope of the parallel code in myReduce(),
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
    T mult = static_cast<T>(multiple);
    for (Iter it = begin; it != end; ++it) {
      MATHICGB_ASSERT(it.index() < colCount());
      // Watch out for overflow here! This is only OK because
      // T is promoting the multiplication to type T.
      mEntries[it.index()] += it.scalar() * mult;
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

    auto begin = matrix.rowBegin(pivotRow);
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

    auto begin = matrixLeft.rowBegin(pivotRow);
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

    auto begin = matrix.rowBegin(pivotRow);
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

    auto begin = matrix.rowBegin(pivotRow);
    SparseMatrix::ColIndex col = begin.index();
    SparseMatrix::Scalar entry = mEntries[col] % modulus;
    mEntries[col] = 0;
    if (entry == 0)
      return;
    SparseMatrix::Scalar reducerLead = begin.scalar();
    ++begin; // can skip first entry as we just set it to zero.
    const auto end = matrix.rowEnd(pivotRow);
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

void myReduce
(SparseMatrix const& toReduceLeft,
 SparseMatrix const& toReduceRight,
 SparseMatrix const& reduceByLeft,
 SparseMatrix const& reduceByRight,
 SparseMatrix::Scalar modulus,
 SparseMatrix& reduced,
 int threadCount) {
  MATHICGB_ASSERT(reduceByLeft.colCount() == reduceByLeft.rowCount());
  const auto pivotCount = reduceByLeft.colCount();
  const auto rowCount = toReduceLeft.rowCount();
  const auto colCountLeft = toReduceLeft.colCount();
  const auto colCountRight = toReduceRight.colCount();

  // ** pre-calculate what rows are pivots for what columns.

  // Store column indexes as the matrix is square anyway (so all indices
  // fit) and we are going to store this as a column index later on.
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

  reduced.clear(colCountRight);

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
	  for (SparseMatrix::ColIndex col = 0; col < colCountRight; ++col) {
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
    dense[row].reset(colCount);
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
          denseRow.normalize(modulus, col);
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
      MATHICGB_ASSERT(dense[row][nextReducers[i].first] == 1); // already normalized
      MATHICGB_ASSERT(reduced.rowCount() == i);
      MATHICGB_ASSERT(!isPivotRow[row]);

      dense[row].appendTo(reduced); // already normalized
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

// doesn't need to be fast.
int integerLog10(size_t val) {
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
    auto it = matrix.rowBegin(row);
    const auto end = matrix.rowEnd(row);
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

void F4MatrixReducer::reduce
(const PolyRing& ring, QuadMatrix& matrix, SparseMatrix& newPivots) {
  MATHICGB_ASSERT(mThreadCount >= 1);
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

  MATHICGB_ASSERT(ring.charac() <= std::numeric_limits<SparseMatrix::Scalar>::max());
  SparseMatrix::Scalar modulus = static_cast<SparseMatrix::Scalar>(ring.charac());

  {

    const SparseMatrix::ColIndex pivotColCount = matrix.topLeft.colCount();

    //SparseMatrix matrixCD;    concatenateMatricesHorizontal    (matrix.bottomLeft, matrix.bottomRight, matrixCD);

    myReduce(matrix.bottomLeft, matrix.bottomRight, matrix.topLeft, matrix.topRight, modulus, newPivots, mThreadCount);
    //    myReduce(matrixCD, matrix.topLeft, matrix.topRight, modulus, newPivots, mThreadCount);    newPivots.trimLeadingZeroColumns(pivotColCount);
  }

  myReduceToEchelonForm5(newPivots, modulus, mThreadCount);
}

F4MatrixReducer::F4MatrixReducer(int threadCount): 
  mThreadCount(std::max(threadCount, 1)) {}
