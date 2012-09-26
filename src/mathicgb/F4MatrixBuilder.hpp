#ifndef F4_MATRIX_BUILDER_GUARD
#define F4_MATRIX_BUILDER_GUARD

#include "QuadMatrixBuilder.hpp"
#include "Poly.hpp"
#include "PolyRing.hpp"
#include "PolyBasis.hpp"
#include <vector>

/** Class for constructing an F4 matrix. This class is reponsible for
  figuring out what matrix to build and then it uses QuadMatrixBuilder
  to create that matrix.

  @todo: this class does not offer exception guarantees. It's just not
  very workable without an RAII monomial handle, so add one of those
  before fixing this. */
class F4MatrixBuilder {
private:
  typedef QuadMatrixBuilder::ColIndex ColIndex;
  typedef QuadMatrixBuilder::LeftRightColIndex LeftRightColIndex;

public:
  F4MatrixBuilder(const PolyBasis& basis):
    mBasis(basis), mBuilder(basis.ring()) {}

  /** Schedules two rows to be added to the matrix whose linear span
    includes the S-polynomial between polyA and polyB. More precisely,
    the two rows represent (B:A)*polyA and (A:B)*polyB where
    A=lead(polyA) and B=lead(polyB). */
  void addTwoRowsForSPairToMatrix(const Poly& polyA, const Poly& polyB) {
    MATHICGB_ASSERT(!polyA.isZero());
    MATHICGB_ASSERT(!polyB.isZero());

    monomial lcm = ring().allocMonomial();
    ring().monomialLeastCommonMultiple
      (polyA.getLeadMonomial(), polyB.getLeadMonomial(), lcm);

    monomial multiple = ring().allocMonomial();

    ring().monomialDivide(polyA.getLeadMonomial(), lcm, multiple);
    addRowToMatrix(multiple, polyA);

    ring().monomialDivide(polyB.getLeadMonomial(), lcm, multiple);
    addRowToMatrix(multiple, polyB);

    ring().freeMonomial(lcm);
    ring().freeMonomial(multiple);
  }

  /** Schedules a row representing multiple*poly to be added to the
      matrix. No ownership is taken, but poly must remain valid until
      the matrix is constructed. multiple is copied so there is no
      requirement there. */
  void addRowToMatrix(const_monomial multiple, const Poly& poly) {
    RowTask task;
    task.useAsReducer = false; // to be updated later
    task.poly = &poly;
    task.multiple = ring().allocMonomial();
    ring().monomialCopy(multiple, task.multiple);
    mTodo.push_back(task);
  }

  /** Builds a matrix and returns it as 4 submatrices A, B, C and
      D. These fit together as
    
      A B
      C D

      where A contains all the pivot rows and columns. The vector tells
      you which monomial each column in B and D represents. This vector
      is sorted in descending order using the order from the basis.

      The matrix contains a reducer/pivot for every monomial that can be
      reduced by the basis and that is present in the matrix. Note that
      the client-added rows also count as reducers so their lead terms
      will not get another reducer added automatically -- specifically,
      adding an S-polynomial will not do what you want because its lead
      term will have no reducer other than itself. Instead, add the two
      polynomials that you would have subtracted from each other to form
      the S-polynomial. */
  void buildMatricesAndClear
  (SparseMatrix& topLeft,
   SparseMatrix& topRight,
   SparseMatrix& bottomLeft,
   SparseMatrix& bottomRight,
   std::vector<monomial> monomialsOfRightColumns)
  {
    // todo: sort input rows by sparsity and/or age.
    // todo: detect and remove duplicate input rows.
    monomial mono = ring().allocMonomial();

    // Decide which input rows are going to be used as reducers and
    // create pivot columns ahead of time for those reducers so that
    // we do not add a reducer for those columns later on.
    typedef std::vector<RowTask>::iterator TaskIter;
    TaskIter end = mTodo.end();
    for (TaskIter it = mTodo.begin(); it != end; ++it) {
      ring().monomialMult(it->multiple, it->poly->getLeadMonomial(), mono);
      LeftRightColIndex leadCol = mBuilder.findColumn(mono);
      it->useAsReducer = !leadCol.valid();
      if (it->useAsReducer) {
        // create column so we know later on that we already have a
        // reducer for this column.
        createOrFindColumnOf(mono);
      }
    }

    // Process pending rows until we are done. Note that the methods
    // we are calling here can add more items to mTodo.
    while (!mTodo.empty()) {
      RowTask task = mTodo.back();
      mTodo.pop_back();
      if (task.useAsReducer)
        appendRowTop(task.multiple, *task.poly);
      else
        appendRowBottom(task.multiple, *task.poly);
    }

    //mBuilder.extractMatricesAndClear
    //(topLeft, topRight, bottomLeft, bottomRight, monomialsOfRightColumns);
  }

  const PolyRing& ring() const {return mBuilder.ring();}

private:
  /** Returns a left or right column that represents mono. Creates a
      new column and schedules a new row to reduce that column if necessary. */
  LeftRightColIndex createOrFindColumnOf(const_monomial mono) {
    LeftRightColIndex colIndex = mBuilder.findColumn(mono);
    if (colIndex.valid())
      return colIndex;

    // mono did not already have a column so look for a reducer
    size_t reducerIndex = mBasis.divisor(mono);
    if (reducerIndex == static_cast<size_t>(-1))
      return LeftRightColIndex(mBuilder.createColumnRight(mono), false);

    // schedule the reducer to be added as a row
    RowTask task;
    task.poly = &mBasis.poly(reducerIndex);
    task.useAsReducer = true;
    task.multiple = ring().allocMonomial();
    ring().monomialDivideToNegative
      (mono, task.poly->getLeadMonomial(), task.multiple);
    mTodo.push_back(task);

    return LeftRightColIndex(mBuilder.createColumnLeft(mono), true);
  }

  void appendRowTop(const_monomial multiple, const Poly& poly) {
    monomial mono = ring().allocMonomial();
    Poly::const_iterator end = poly.end();
    for (Poly::const_iterator it = poly.begin(); it != end; ++it) {
      ring().monomialMult(it.getMonomial(), multiple, mono);
      mBuilder.appendEntryTop(createOrFindColumnOf(mono), it.getCoefficient());
    }
    ring().freeMonomial(mono);
  }

  void appendRowBottom(const_monomial multiple, const Poly& poly) {
    monomial mono = ring().allocMonomial();
    Poly::const_iterator end = poly.end();
    for (Poly::const_iterator it = poly.begin(); it != end; ++it) {
      ring().monomialMult(it.getMonomial(), multiple, mono);
      mBuilder.appendEntryBottom
        (createOrFindColumnOf(mono), it.getCoefficient());
    }
    ring().freeMonomial(mono);
  }

  /// Represents the task of adding a row representing poly*multiple.
  struct RowTask {
    bool useAsReducer; // if true: put in top part of matrix
    const Poly* poly;
    monomial multiple;
  };
  std::vector<RowTask> mTodo;
  const PolyBasis& mBasis;
  QuadMatrixBuilder mBuilder;
};

#endif
