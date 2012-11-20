#ifndef MATHICGB_SPARSE_MATRIX_GUARD
#define MATHICGB_SPARSE_MATRIX_GUARD

#include "RawVector.hpp"
#include "PolyRing.hpp"
#include <mathic.h>
#include <vector>
#include <ostream>
#include <limits>
#include <sstream>
#include <string>
#include <iterator>
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
*/
class SparseMatrix {
public:
  typedef size_t RowIndex;
  typedef uint32 ColIndex;
  typedef uint16 Scalar;
  class ConstRowIterator;

  /// Construct a matrix with no rows.
  SparseMatrix(const size_t memoryQuantum = 0):
    mMemoryQuantum(memoryQuantum)
  {}

  SparseMatrix(SparseMatrix&& matrix):
    mRows(std::move(matrix.mRows)),
    mBlock(std::move(matrix.mBlock)),
    mMemoryQuantum(matrix.mMemoryQuantum)
  {
  }

  SparseMatrix& operator=(SparseMatrix&& matrix) {
    this->~SparseMatrix();
    new (this) SparseMatrix(std::move(matrix));
    return *this;
  }

  SparseMatrix(const SparseMatrix& matrix) {
    *this = matrix;
  }

  ~SparseMatrix() {clear();}

  SparseMatrix& operator=(const SparseMatrix&);
  void swap(SparseMatrix& matrix);

  bool operator==(const SparseMatrix& matrix) const;
  bool operator!=(const SparseMatrix& matrix) const {
    return !(*this == matrix);
  }

  // Removes all rows from *this.
  void clear();

  /// Appends the rows from matrix to this object. Avoids most of the copies
  /// that would otherwise be required for a big matrix insert by taking
  /// the memory out of matrix.
  void takeRowsFrom(SparseMatrix&& matrix);

  RowIndex rowCount() const {return mRows.size();}
  ColIndex computeColCount() const;
  size_t memoryQuantum() const {return mMemoryQuantum;}

  /// Returns the number of entries in the whole matrix. Is not constant time
  /// so avoid calling too many times.
  size_t entryCount() const {
    size_t count = 0;
    const Block* block = &mBlock;
    for (; block != 0; block = block->mPreviousBlock)
      count += block->mColIndices.size();
    return count;
  }

  /// Returns the number of bytes of memory allocated by this object. Is not
  /// constant time so avoid calling too many times.
  size_t memoryUse() const;
  size_t memoryUseTrimmed() const;

   /// Returns the number of entries in the given row.
  ColIndex entryCountInRow(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    return mRows[row].size();
  }

  /// Returns true if the given row has no entries.
  bool emptyRow(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    return mRows[row].empty();
  }

  ConstRowIterator rowBegin(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    const Row& r = mRows[row];
    return ConstRowIterator(r.mIndicesBegin, r.mScalarsBegin);
  }

  ConstRowIterator rowEnd(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    const Row& r = mRows[row];
    return ConstRowIterator(r.mIndicesEnd, r.mScalarsEnd);
  }

  /// Returns the index of the first entry in the given row. This is
  /// the first entry that you added to the row - so not necessarily the
  /// minimum column index in that row. The row in question must have at
  /// least one entry.
  ColIndex leadCol(RowIndex row) const {
    MATHICGB_ASSERT(row < rowCount());
    MATHICGB_ASSERT(!emptyRow(row));
    return *mRows[row].mIndicesBegin;
  }

  /// Prints the matrix in a human readable format to out.
  void print(std::ostream& out) const;

  std::string toString() const;



  /// Removes the leading trimThisMany columns. The columns are
  /// removed by replacing all column indices col by col -
  /// trimThisMany. No entry can have a column index less than
  /// trimThisMany, even if the scalar of that entry is set to zero.
  void trimLeadingZeroColumns(ColIndex trimThisMany);

  /// Ensure that there is enough space for at least freeCount additional
  /// entries without needing to allocate more memory for entries.
  /// Pending entries that are not fixed into a row yet do not count as
  /// free for this calculation.
  void reserveFreeEntries(size_t freeCount);

  /// Preallocate space for at least count rows.
  void reserveRows(size_t count) {mRows.reserve(count);}

  /// Adds a new row that contains all terms that have been appended
  /// since the last time a row was added or the matrix was created.
  void rowDone() {
    MATHICGB_ASSERT(mBlock.mColIndices.size() == mBlock.mScalars.size());
    Row row;
    row.mIndicesEnd = mBlock.mColIndices.end();
    row.mScalarsEnd = mBlock.mScalars.end();
    if (mBlock.mHasNoRows) {
      row.mIndicesBegin = mBlock.mColIndices.begin();
      row.mScalarsBegin = mBlock.mScalars.begin();
      mBlock.mHasNoRows = false;
    } else {
      row.mIndicesBegin = mRows.back().mIndicesEnd;
      row.mScalarsBegin = mRows.back().mScalarsEnd;
    }
    mRows.push_back(row);
  }

  /// Appends an entry to the matrix. Will not appear in the matrix
  /// until rowDone is called. Do not call other methods that add rows
  /// after calling this method until rowDone has been called.
  inline void appendEntry(ColIndex colIndex, Scalar scalar) {
    MATHICGB_ASSERT(mBlock.mColIndices.size() == mBlock.mScalars.size());

    MATHICGB_ASSERT(mBlock.mScalars.atCapacity() ==
      mBlock.mColIndices.atCapacity());
    if (mBlock.mScalars.atCapacity())
      growEntryCapacity();
    MATHICGB_ASSERT(!mBlock.mScalars.atCapacity());
    MATHICGB_ASSERT(!mBlock.mColIndices.atCapacity());

    mBlock.mColIndices.rawPushBack(colIndex);
    mBlock.mScalars.rawPushBack(scalar);

    MATHICGB_ASSERT(mBlock.mColIndices.size() == mBlock.mScalars.size());
  }

  void appendRowAndNormalize(const SparseMatrix& matrix, RowIndex row, Scalar modulus);
  
  void appendRow(const SparseMatrix& matrix, RowIndex row);

  void appendRowWithModulus(const std::vector<uint64>& v, Scalar modulus);
  
  void appendRow(const std::vector<uint64>& v, ColIndex leadCol = 0);

  void appendRowWithModulusNormalized(const std::vector<uint64>& v, Scalar modulus);

  // Returns true if the row was non-zero. Otherwise the row was not
  // appended.
  bool appendRowWithModulusIfNonZero(const std::vector<uint64>& v, Scalar modulus);

  /// Replaces all column indices i with colMap[i].
  void applyColumnMap(const std::vector<ColIndex>& colMap);

  /// Let poly be the dot product of colMonomials and the given row.
  void rowToPolynomial(
    RowIndex row,
    const std::vector<monomial>& colMonomials,
    Poly& poly);

  /// Reorders the rows so that the index of the leading column in
  /// each row is weakly increasing going from top to bottom. Quite
  /// slow and it makes a copy internally.
  void sortRowsByIncreasingPivots();

  // Write *this and modulus to file.
  void write(Scalar modulus, FILE* file) const;

  // Set *this to a matrix read from file and return the modulus from the file.
  Scalar read(FILE* file);

  /// Iterates through the entries in a row.
  class ConstRowIterator {
  public:
    typedef const std::pair<ColIndex, Scalar> value_type;
	typedef ptrdiff_t difference_type;
    typedef size_t distance_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::random_access_iterator_tag iterator_category;

    ConstRowIterator& operator++() {
      ++mScalarIt;
      ++mColIndexIt;
      return *this;
    }

    value_type operator*() const {return value_type(index(), scalar());}

    bool operator<(const ConstRowIterator& it) const {
      return mColIndexIt < it.mColIndexIt;
    }

    difference_type operator-(const ConstRowIterator& it) const {
      return mColIndexIt - it.mColIndexIt;
    }

    bool operator==(const ConstRowIterator& it) const {
      return mColIndexIt == it.mColIndexIt;
    }

    bool operator!=(const ConstRowIterator& it) const {return !(*this == it);}
    const Scalar& scalar() const {return *mScalarIt;}
    const ColIndex& index() const {return *mColIndexIt;}

  private:
    friend class SparseMatrix;
    ConstRowIterator(
      const ColIndex* const indicesIt,
      const Scalar* const scalarIt
    ):
      mColIndexIt(indicesIt),
      mScalarIt(scalarIt)
    {
    }

    const ColIndex* mColIndexIt;
    const Scalar* mScalarIt;
  };

private:
  MATHICGB_NO_INLINE void growEntryCapacity();

  /// Contains information about a row in the matrix.
  struct Row {
    Row(): mScalarsBegin(0), mScalarsEnd(0), mIndicesBegin(0), mIndicesEnd(0) {}

    Scalar* mScalarsBegin;
    Scalar* mScalarsEnd;
    ColIndex* mIndicesBegin;
    ColIndex* mIndicesEnd;

    bool empty() const {return mIndicesBegin == mIndicesEnd;}
    ColIndex size() const {
      return static_cast<ColIndex>(std::distance(mIndicesBegin, mIndicesEnd));
    }
  };
  std::vector<Row> mRows;

  /// Memory is allocated a block at a time. This avoids the need for copying
  /// that a std::vector normally does on reallocation. Believe it or not,
  /// copying sparse matrix memory due to reallocation was accounting for 5%
  /// of the running time before this change.
  struct Block {
    Block(): mPreviousBlock(0), mHasNoRows(true) {}
    Block(Block&& block):
      mColIndices(std::move(block.mColIndices)),
      mScalars(std::move(block.mScalars)),
      mPreviousBlock(block.mPreviousBlock),
      mHasNoRows(block.mHasNoRows) 
    {
      block.mPreviousBlock = 0;
      block.mHasNoRows = true;
    }

    void swap(Block& block) {
      std::swap(mColIndices, block.mColIndices);
      std::swap(mScalars, block.mScalars);
      std::swap(mPreviousBlock, block.mPreviousBlock);
      std::swap(mHasNoRows, block.mHasNoRows);
    }

    Block& operator=(Block&& block) {
      this->~Block();
      new (this) Block(std::move(block));
      return *this;
    }

    size_t memoryUse() const;
    size_t memoryUseTrimmed() const;

    /// We need a RawVector here to tie the checks for the need to reallocate
    /// together between mColIndices and mEntries. We only need to check
    /// the capacity once, which, believe it or not, is a significant performance
    /// win. Not least because it decreases the amount of code and therefore
    /// causes better compiler inlining decisions.
    RawVector<ColIndex> mColIndices;
    RawVector<Scalar> mScalars;
    Block* mPreviousBlock; /// is null if there are no previous blocks
    bool mHasNoRows; /// true if no rows have been made from this block yet
  };
  Block mBlock;
  size_t mMemoryQuantum;
};

inline void swap(SparseMatrix& a, SparseMatrix& b) {
  a.swap(b);
}

std::ostream& operator<<(std::ostream& out, const SparseMatrix& matrix);

#endif
