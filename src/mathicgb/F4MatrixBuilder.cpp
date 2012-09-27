#include "stdinc.h"
#include "F4MatrixBuilder.hpp"

F4MatrixBuilder::F4MatrixBuilder(const PolyBasis& basis):
  mBasis(basis), mBuilder(basis.ring()) {}

void F4MatrixBuilder::addTwoRowsForSPairToMatrix(const Poly& polyA, const Poly& polyB) {
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

void F4MatrixBuilder::addRowToMatrix
(const_monomial multiple, const Poly& poly) {
  RowTask task;
  task.useAsReducer = false; // to be updated later
  task.poly = &poly;

  task.multiple = ring().allocMonomial();
  ring().monomialCopy(multiple, task.multiple);
  mTodo.push_back(task);
}

void F4MatrixBuilder::buildMatrixAndClear(QuadMatrix& matrix) {
  // todo: detect and remove duplicate input rows.
  // todo: prefer sparse/old reducers among the inputs.
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
      mBuilder.createColumnLeft(mono);
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

  mBuilder.buildMatrixAndClear(matrix);
}

F4MatrixBuilder::LeftRightColIndex
F4MatrixBuilder::createOrFindColumnOf(const_monomial mono) {
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

void F4MatrixBuilder::appendRowTop(const_monomial multiple, const Poly& poly) {
  monomial mono = ring().allocMonomial();
  Poly::const_iterator end = poly.end();
  for (Poly::const_iterator it = poly.begin(); it != end; ++it) {
    ring().monomialMult(it.getMonomial(), multiple, mono);
    mBuilder.appendEntryTop(createOrFindColumnOf(mono), it.getCoefficient());
  }
  mBuilder.rowDoneTopLeftAndRight();
  ring().freeMonomial(mono);
}

void F4MatrixBuilder::appendRowBottom
(const_monomial multiple, const Poly& poly) {
  monomial mono = ring().allocMonomial();
  Poly::const_iterator end = poly.end();
  for (Poly::const_iterator it = poly.begin(); it != end; ++it) {
    ring().monomialMult(it.getMonomial(), multiple, mono);
    mBuilder.appendEntryBottom
      (createOrFindColumnOf(mono), it.getCoefficient());
  }
  mBuilder.rowDoneBottomLeftAndRight();
  ring().freeMonomial(mono);
}
