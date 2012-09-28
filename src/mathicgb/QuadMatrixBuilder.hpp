#ifndef MATHICGB_QUAD_MATRIX_BUILDER_GUARD
#define MATHICGB_QUAD_MATRIX_BUILDER_GUARD

#define MATHICGB_USE_QUADMATRIX_STD_HASH

#include "SparseMatrix.hpp"
#include "PolyRing.hpp"
#include <vector>
#include <map>
#include <string>
#include <ostream>
#ifdef MATHICGB_USE_QUADMATRIX_STD_HASH
#include <unordered_map>
#endif
class FreeModuleOrder;
class QuadMatrix;

/** Builder for QuadMatrix. This is not quite the builder pattern in
  that the interface is not virtual and the implementation cannot be
  swapped out - it only follows the builder pattern in that it is a
  class that allows step-wise construction of a final product. */
class QuadMatrixBuilder {
 public:
  typedef SparseMatrix::RowIndex RowIndex;
  typedef SparseMatrix::ColIndex ColIndex;
  typedef SparseMatrix::Scalar Scalar;

  QuadMatrixBuilder(const PolyRing& ring):
#ifndef MATHICGB_USE_QUADMATRIX_STD_HASH
    mMonomialToCol(ArbitraryOrdering(ring)) {}
#else
  mMonomialToCol(100, Hash(ring), Equal(ring)) {}
#endif

  /// The index of a column that can be either on the left or the
  /// right side. The largest representable ColIndex is an invalid
  /// index. This is the default value. The only allowed method to
  /// call for an invalid index is valid().
  class LeftRightColIndex {
  public:
    LeftRightColIndex():
      mRawIndex(std::numeric_limits<ColIndex>::max()), mLeft(false) {}
    LeftRightColIndex(ColIndex index, bool left):
      mRawIndex(index), mLeft(left) {
    }

    ColIndex leftIndex() const {
      MATHICGB_ASSERT(left());
      return index();
    }

    ColIndex rightIndex() const {
      MATHICGB_ASSERT(right());
      return index();
    }

    /// Use leftIndex() or rightIndex() instead if you know what side
    /// you are expecting, as check the side in assert mode.
    ColIndex index() const {
      MATHICGB_ASSERT(valid());
      return mRawIndex;
    }

    bool left() const {
      MATHICGB_ASSERT(valid());
      return mLeft;
    }

    bool right() const {
      MATHICGB_ASSERT(valid());
      return !left();
    }

    bool valid() const {
      return mRawIndex != std::numeric_limits<ColIndex>::max();
    }

  private:
    ColIndex mRawIndex;
    bool mLeft;
  };

  // **** Appending entries to top matrices.
  // Same interface as SparseMatrix except with two matrices and here
  // you have to create columns before you can use them.

  void appendEntryTopLeft(ColIndex col, Scalar scalar) {
    MATHICGB_ASSERT(col < leftColCount());
    mTopLeft.appendEntry(col, scalar);
  }

  void appendEntryTopRight(ColIndex col, Scalar scalar) {
    MATHICGB_ASSERT(col < rightColCount());
    mTopRight.appendEntry(col, scalar);
  }

  void appendEntryTop(LeftRightColIndex col, Scalar scalar) {
    MATHICGB_ASSERT(col.valid());
    if (col.left())
      appendEntryTopLeft(col.leftIndex(), scalar);
    else
      appendEntryTopRight(col.rightIndex(), scalar);
  }

  void rowDoneTopLeftAndRight() {
    mTopLeft.rowDone();
    mTopRight.rowDone();
  }

  // **** Appending entries to bottom matrices
  // Same interface as SparseMatrix except with two matrices and here
  // you have to create columns before you can use them.

  void appendEntryBottomLeft(ColIndex col, Scalar scalar) {
    MATHICGB_ASSERT(col < leftColCount());
    mBottomLeft.appendEntry(col, scalar);
  }

  void appendEntryBottomRight(ColIndex col, Scalar scalar) {
    MATHICGB_ASSERT(col < rightColCount());
    mBottomRight.appendEntry(col, scalar);
  }

  void appendEntryBottom(LeftRightColIndex col, Scalar scalar) {
    MATHICGB_ASSERT(col.valid());
    if (col.left())
      appendEntryBottomLeft(col.leftIndex(), scalar);
    else
      appendEntryBottomRight(col.rightIndex(), scalar);
  }

  void rowDoneBottomLeftAndRight() {
    mBottomLeft.rowDone();
    mBottomRight.rowDone();
  }

  // *** Creating and reordering columns
  // Unlike the interface for SparseMatrix, here you have to create
  // columns before you can append entries in those columns. All
  // passed in monomials are copied so that ownership of the memory is
  // not taken over.

  /** Creates a new column associated to the monomial
    monomialToBeCopied to the left matrices. There must not already
    exist a column for this monomial on the left or on the right. */
  ColIndex createColumnLeft(const_monomial monomialToBeCopied);

  /** Creates a new column associated to the monomial monomialToBeCopied
    to the right matrices. There must not already exist a column for
    this monomial on the left or on the right. */
  ColIndex createColumnRight(const_monomial monomialToBeCopied);

  /** Sorts the left columns to be decreasing with respect to
    order. Also updates the column indices already in the matrix to
    reflect the new ordering. */
  void sortColumnsLeft(const FreeModuleOrder& order);

  /** Sorts the right columns to be decreasing with respect to
    order. Also updates the column indices already in the matrix to
    reflect the new ordering. */
  void sortColumnsRight(const FreeModuleOrder& order);


  // *** Querying columns

  // Returns a column for the findThis monomial. Searches on both the
  // left and right side. Returns an invalid index if no such column
  // exists.
  LeftRightColIndex findColumn(const_monomial findThis) const {
    MonomialToColType::const_iterator it = mMonomialToCol.find(findThis);
    if (it != mMonomialToCol.end())
      return it->second;
    else
      return LeftRightColIndex();
  }

  const_monomial monomialOfLeftCol(ColIndex col) const {
    MATHICGB_ASSERT(col < mMonomialsLeft.size());
    return mMonomialsLeft[col];
  }

  const_monomial monomialOfRightCol(ColIndex col) const {
    MATHICGB_ASSERT(col < mMonomialsRight.size());
    return mMonomialsRight[col];
  }

  const_monomial monomialOfCol(LeftRightColIndex col) const {
    MATHICGB_ASSERT(col.valid());
    if (col.left())
      return monomialOfLeftCol(col.leftIndex());
    else
      return monomialOfRightCol(col.rightIndex());
  }

  const SparseMatrix& topLeft() const {return mTopLeft;}
  const SparseMatrix& topRight() const {return mTopRight;}
  const SparseMatrix& bottomLeft() const {return mBottomLeft;}
  const SparseMatrix& bottomRight() const {return mBottomRight;}

  // String representation intended for debugging.
  void print(std::ostream& out) const;

  // String representation intended for debugging.
  std::string toString() const;

  const PolyRing& ring() const {
    // The key comparer object already has a ring reference - we might
    // as well use that one instead of adding another reference to
    // this object.
#ifndef MATHICGB_USE_QUADMATRIX_STD_HASH
    return mMonomialToCol.key_comp().ring();
#else
    return mMonomialToCol.key_eq().ring();
#endif
  }

  ColIndex leftColCount() const {
    MATHICGB_ASSERT(topLeft().colCount() == bottomLeft().colCount());
    return topLeft().colCount();
  }

  ColIndex rightColCount() const {
    MATHICGB_ASSERT(topRight().colCount() == bottomRight().colCount());
    return topRight().colCount();
  }

  /// Puts the built matrix into out and sets the builder to a state
  /// with no columns and no rows.
  void buildMatrixAndClear(QuadMatrix& out);

private:
  // Store the monomials for each left and right column respectivewy.
  typedef std::vector<monomial> MonomialsType;
  MonomialsType mMonomialsLeft;
  MonomialsType mMonomialsRight;

#ifndef MATHICGB_USE_QUADMATRIX_STD_HASH
  /// We need SOME ordering to make std::map work.
  class ArbitraryOrdering {
  public:
    ArbitraryOrdering(const PolyRing& ring): mRing(ring) {}
    bool operator()(const_monomial a, const_monomial b) const {
      return mRing.monomialLT(a, b);
    }
    const PolyRing& ring() const {return mRing;}

  private:
    const PolyRing& mRing;
  };
  /// Allows fast determination of which column has a given monomial.
  typedef std::map<const_monomial, LeftRightColIndex, ArbitraryOrdering>
    MonomialToColType;
  MonomialToColType mMonomialToCol;
#else
  struct Hash {
  public:
    Hash(const PolyRing& ring): mRing(ring) {}
    size_t operator()(const_monomial m) const {
      MATHICGB_ASSERT(mRing.hashValid(m));
      return mRing.monomialHashValue(m);
    }
    const PolyRing& ring() const {return mRing;}

  private:
    const PolyRing& mRing;
  };
  struct Equal {
  public:
    Equal(const PolyRing& ring): mRing(ring) {}
    size_t operator()(const_monomial a, const_monomial b) const {
      return mRing.monomialEQ(a, b);
    }
    const PolyRing& ring() const {return mRing;}

  private:
    const PolyRing& mRing;
  };

  typedef std::unordered_map<const_monomial, LeftRightColIndex, Hash, Equal> 
    MonomialToColType;
  MonomialToColType mMonomialToCol;
#endif

  SparseMatrix mTopLeft;
  SparseMatrix mTopRight;
  SparseMatrix mBottomLeft;
  SparseMatrix mBottomRight;
};

#endif
