#ifndef _bit_triangle_h_
#define _bit_triangle_h_

#include <vector>

// Object that stores a triangular 2-dimensional array of bits. For example:
//
// row
//  3|      
//  2|      1
//  1|    0 1
//  0|  0 0 0
//    -------
//    0 1 2 3 column
//
// A bit is addressed by a pair (column, row) where the column goes first.
// All valid address pairs have 0 <= row < column < columnCount().
class BitTriangle {
public:
  // Returns how many columns the triangle has
  size_t columnCount() const {return mColumns.size();}

  // Returns true if there are no columns in the triangle
  bool empty() const {return mColumns.empty();}

  // Adds a new column of the triangle. This increases columnCount() by
  // one, and the index of the new column is the previous value of
  // columnCount(). The new bits are all set to false initially.
  void addColumn() {
    size_t oldSize = mColumns.size();
    mColumns.resize(oldSize + 1);
    mColumns[oldSize].resize(oldSize);
  }

  // Returns the bit in the given column and row. As this is a triangle it
  // must be true that column >= row.
  bool bit(size_t column, size_t row) const {
    ASSERT(column < columnCount());
    ASSERT(row < column);
    return mColumns[column][row];
  }

  // As bit, but uses max(x,y) as the column and min(x,y) as the row.
  bool bitUnordered(size_t x, size_t y) const {
    ASSERT(x < columnCount());
    ASSERT(y < columnCount());
    ASSERT(x != y);
    if (x < y)
      std::swap(x, y);
    return bit(x, y);
  }

  // Returns the bit in the given column and row. As this is a triangle it
  // must be true that column >= row.
  void setBit(size_t column, size_t row, bool value) {
    ASSERT(column < columnCount());
    ASSERT(row < column);
    mColumns[column][row] = value;
  }

  // As setBit, but uses max(x,y) as the column and min(x,y) as the row.
  void setBitUnordered(size_t x, size_t y, bool value) {
    ASSERT(x < columnCount());
    ASSERT(y < columnCount());
    ASSERT(x != y);
    if (x < y)
      std::swap(x, y);
    setBit(x, y, value);
  }

  size_t getMemoryUse() const;

private:
  std::vector<std::vector<bool> > mColumns;
};

#endif

