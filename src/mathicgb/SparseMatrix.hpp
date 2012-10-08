#ifndef MATHICGB_SPARSE_MATRIX_GUARD
#define MATHICGB_SPARSE_MATRIX_GUARD

#include "PolyRing.hpp"
#include <vector>
#include <ostream>
#include <limits>
#include <sstream>
#include <string>
class Poly;

/** A class that implements a sparse matrix.

These are the mathematical concepts involved:
  Sparse matrix: a sequence of sparse rows.
  Sparse row: a seqence of entries.
  Entry: a pair (i,s) where i is a column index and s is a scalar.

You add a row by adding all entries in the row and then calling
rowDone(). You cannot add entries to a row once it has been
created, so in that sense this class is append-only. However, you
are free to change the indices and the scalars in the entries that
are already there. Entries are not automatically reordered by this
class, so your rows will be in increasing order of index only if
you make them like that.

Adding an entry or a row can invalidate all pointers/references to
entries in the matrix and all iterators. This is true even if the
entry has been added but it has not been put in a new row yet by
calling rowDone.

There is no special treatment of entries whose scalar is zero. For
example they still count as entries in relation to entryCount().

Currently this is not a template class so you can get by without
using the typedefs offered, for example using uint16 instead of
SparseMatrix::Scalar. Please use the typedefs to make it easier to
support a wider range of types of matrices in future.

There is no need to specify the number of columns ahead of time. Any
column index within the range of the ColIndex type can be used. The
SparseMatrix automatically keeps track of the number of columns in the
matrix.
*/
class SparseMatrix {
 public:
  typedef size_t RowIndex;
  typedef uint32 ColIndex;
  typedef uint16 Scalar;

  /** Preallocate space for at least count entries. */
  void reserveEntries(size_t count) {
    mEntries.reserve(count);
    mColIndices.reserve(count);
  }

  /** Preallocate space for at least count rows. */
  void reserveRows(size_t count) {
    mRowOffsets.reserve(count + 1);
  }

  /** Returns the index of the first entry in the given row. This is
    the first entry that you added to the row - so not necessarily the
    minimum column index in that row. The row in question must have at
    least one entry. */
  ColIndex leadCol(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    MATHICGB_ASSERT(!emptyRow(row));
    return mColIndices[mRowOffsets[row]];
  }

  /** Returns true if the given row has no entries. */
  bool emptyRow(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    return mRowOffsets[row] == mRowOffsets[row + 1];
  }

  /** Removes the leading trimThisMany columns. The columns are
    removed by replacing all column indices col by col -
    trimThisMany. No entry can have a column index less than
    trimThisMany, even if the scalar of that entry is set to zero. */
  void trimLeadingZeroColumns(ColIndex trimThisMany) {
    MATHICGB_ASSERT(trimThisMany <= colCount());
    typedef std::vector<ColIndex>::iterator Iter;
    Iter end = mColIndices.end();
    Iter it = mColIndices.begin();
    for (; it != end; ++it) {
      MATHICGB_ASSERT(*it >= trimThisMany);
      *it -= trimThisMany;
    }
    mColCount -= trimThisMany;
  }

  /** Construct a new zero by zero sparse matrix. */
  SparseMatrix(): mColCount(0) {
    mRowOffsets.push_back(0);
  }

  /** Returns the number of entries in the given row. */
  size_t entryCountInRow(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    return mRowOffsets[row + 1] - mRowOffsets[row];
  }

  /** Returns the number of entries in the whole matrix. */
  size_t entryCount() const {
    return mEntries.size();
  }

  /** Prints the matrix in a human readable format to out. */
  void print(std::ostream& out) const {
    if (rowCount() == 0)
      out << "matrix with no rows\n";
    for (RowIndex row = 0; row < rowCount(); ++row) {
      out << row << ':';
      RowIterator end = rowEnd(row);
      for (RowIterator it = rowBegin(row); it != end; ++it) {
        MATHICGB_ASSERT(it.index() < colCount());
        out << ' ' << it.index() << '#' << it.scalar();
      }
      out << '\n';
    }
  }

  std::string toString() const {
    std::ostringstream out;
    print(out);
    return out.str();
  }


  /** Adds a new row that contains all terms that have been appended
    since the last time a row was added or the matrix was created. */
  void rowDone() {
    MATHICGB_ASSERT(mColIndices.size() == mEntries.size());

    mRowOffsets.push_back(mColIndices.size());
  }

  /** Appends an entry to the matrix. Will not appear in the matrix
    until rowDone is called. Do not call other methods that add rows
    after calling this method until rowDone has been called. */
  void appendEntry(ColIndex colIndex, Scalar scalar) {
    MATHICGB_ASSERT(mColIndices.size() == mEntries.size());
    MATHICGB_ASSERT(colIndex < std::numeric_limits<ColIndex>::max());

    mColIndices.push_back(colIndex);
    mEntries.push_back(scalar);
    if (colIndex + 1 > mColCount)
      mColCount = colIndex + 1;

    MATHICGB_ASSERT(mColIndices.size() == mEntries.size());
  }

  
  void appendRowAndNormalize(const SparseMatrix& matrix, RowIndex row, Scalar modulus) {
    MATHICGB_ASSERT(row < matrix.rowCount());
    RowIterator it = matrix.rowBegin(row);
    RowIterator end = matrix.rowEnd(row);
    if (it != end) {
      appendEntry(it.index(), 1);
      Scalar lead = it.scalar();
      ++it;
      if (it != end) {
        Scalar inverse = modularInverse(lead, modulus);
        do {
          uint32 prod = static_cast<uint32>(inverse) * it.scalar();
          uint16 prodMod = static_cast<uint16>(prod % modulus);
          appendEntry(it.index(), prodMod);
          ++it;
        } while (it != end);
      }
    }
    rowDone();
  }
  
  void appendRow(const SparseMatrix& matrix, RowIndex row) {
    MATHICGB_ASSERT(row < matrix.rowCount());
    RowIterator it = matrix.rowBegin(row);
    RowIterator end = matrix.rowEnd(row);
    for (; it != end; ++it)
      appendEntry(it.index(), it.scalar());
    rowDone();
  }

  void setToAndTakeMemory(std::vector<ColIndex>& indices,
                          std::vector<Scalar>& entries,
                          std::vector<ColIndex>& sizes,
                          ColIndex colCount) {
    MATHICGB_ASSERT(indices.size() == entries.size());

    mColCount = colCount;

    // deallocate old memory
    std::vector<ColIndex>().swap(mColIndices);
    std::vector<Scalar>().swap(mEntries);
    std::vector<size_t>().swap(mRowOffsets);

    // take new memory
    mColIndices.swap(indices);
    mEntries.swap(entries);

    // compute row offsets
    const size_t newRowCount = sizes.size();
    mRowOffsets.reserve(newRowCount + 1);
    mRowOffsets.push_back(0);
    for (size_t row = 0; row < newRowCount; ++row)
      mRowOffsets.push_back(mRowOffsets.back() + sizes[row]);

    MATHICGB_ASSERT(mRowOffsets.size() == sizes.size() + 1);
    MATHICGB_ASSERT(mRowOffsets.front() == 0);
    MATHICGB_ASSERT(mRowOffsets.back() == mColIndices.size());

    std::vector<ColIndex>().swap(sizes); // deallocate old memory
  }
  
  void swap(SparseMatrix& matrix) {
    mColIndices.swap(matrix.mColIndices);
    mEntries.swap(matrix.mEntries);
    mRowOffsets.swap(matrix.mRowOffsets);
    std::swap(mColCount, matrix.mColCount);
  }

  void clear(ColIndex newColCount = 0) {
    mColIndices.clear();
    mEntries.clear();
    mRowOffsets.clear();
    mRowOffsets.push_back(0);
    mColCount = newColCount;
  }
  
  class RowIterator {
  public:
    RowIterator& operator++() {
      ++mOffset;
      return *this;
    }
    bool operator==(const RowIterator& it) const {return mOffset == it.mOffset;}
    bool operator!=(const RowIterator& it) const {return !(*this == it);}
    Scalar scalar() {return mMatrix.scalarAtOffset(mOffset);}
    ColIndex index() {return mMatrix.indexAtOffset(mOffset);}
    
  private:
    friend class SparseMatrix;
    RowIterator(const SparseMatrix& matrix, size_t offset): mMatrix(matrix), mOffset(offset) {}
 
    const SparseMatrix& mMatrix;
    size_t mOffset;      
  };
  
  RowIndex rowCount() const {
    MATHICGB_ASSERT(!mRowOffsets.empty());
    return mRowOffsets.size() - 1;
  }

  ColIndex colCount() const {
    return mColCount;
  }
  
  RowIterator rowBegin(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    return RowIterator(*this, mRowOffsets[row]);
  }

  RowIterator rowEnd(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    return RowIterator(*this, mRowOffsets[row + 1]);
  }
  
  void ensureAtLeastThisManyColumns(ColIndex count) {
    if (count > colCount())
      mColCount = count;
  }

  void appendRowWithModulus(std::vector<uint64> const& v, Scalar modulus) {
    MATHICGB_ASSERT(v.size() == colCount());
    ColIndex count = colCount();
    for (ColIndex col = 0; col < count; ++col) {
      Scalar scalar = static_cast<Scalar>(v[col] % modulus);
      if (scalar != 0)
        appendEntry(col, scalar);
    }
    rowDone();
  }
  
  void appendRow(std::vector<uint64> const& v, ColIndex leadCol = 0) {
    MATHICGB_ASSERT(v.size() == colCount());
#ifdef MATHICGB_DEBUG
    for (ColIndex col = leadCol; col < leadCol; ++col) {
      MATHICGB_ASSERT(v[col] == 0);
    }
#endif

    ColIndex count = colCount();
    for (ColIndex col = leadCol; col < count; ++col) {
	  MATHICGB_ASSERT(v[col] < std::numeric_limits<Scalar>::max());
      if (v[col] != 0)
        appendEntry(col, static_cast<Scalar>(v[col]));
	}
    rowDone();
  }

  void appendRowWithModulusNormalized(std::vector<uint64> const& v, Scalar modulus) {
    MATHICGB_ASSERT(v.size() == colCount());
    ColIndex count = colCount();
    uint16 multiply = 1;
   
    bool first = true;
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

  // Returns true if the row was non-zero. Otherwise the row was not
  // appended.
  bool appendRowWithModulusIfNonZero(std::vector<uint64> const& v, Scalar modulus) {
    appendRowWithModulus(v, modulus);
    MATHICGB_ASSERT(mRowOffsets.size() >= 2);
    std::vector<size_t>::const_iterator it = mRowOffsets.end();
    --it;
    size_t last = *it;
    --it;
    if (last == *it) {
      mRowOffsets.pop_back();
      return false;
    } else
      return true;
  }

  bool operator==(const SparseMatrix& mat) const {
    return mColCount == mat.mColCount &&
      mColIndices == mat.mColIndices &&
      mEntries == mat.mEntries &&
      mRowOffsets == mat.mRowOffsets;
  }

  bool operator!=(const SparseMatrix& mat) const {
    return !(*this == mat);
  }

  const std::vector<Scalar>& entries() const {
    return mEntries;
  }

  const std::vector<ColIndex>& colIndices() const {
    return mColIndices;
  }
  
  typedef std::vector<ColIndex>::iterator AllColIndicesIterator;
  AllColIndicesIterator allColIndicesBegin() {return mColIndices.begin();}
  AllColIndicesIterator allColIndicesEnd() {return mColIndices.end();}

  /// Let poly be the dot product of colMonomials and the given row.
  void rowToPolynomial
  (RowIndex row, std::vector<monomial> colMonomials, Poly& poly);

  /// Reorders the rows so that the index of the leading column in
  /// each row is weakly increasing going from top to bottom. Quite
  /// slow and it makes a copy internally.
  void sortRowsByIncreasingPivots();
  
private:
  friend class RowIterator;

  ColIndex indexAtOffset(size_t offset) const {
    MATHICGB_ASSERT(offset < mColIndices.size());
    return mColIndices[offset];
  }

  Scalar scalarAtOffset(size_t offset) const {
    MATHICGB_ASSERT(offset < mEntries.size());
    return mEntries[offset];
  }

  std::vector<ColIndex> mColIndices;
  std::vector<Scalar> mEntries;
  std::vector<size_t> mRowOffsets;  
  ColIndex mColCount;
};

std::ostream& operator<<(std::ostream& out, const SparseMatrix& matrix);

#endif
