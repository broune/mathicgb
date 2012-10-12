#include "stdinc.h"
#include "F4MatrixBuilder.hpp"

F4MatrixBuilder::F4MatrixBuilder(const PolyBasis& basis):
  mBasis(basis), mBuilder(basis.ring()) {}

void F4MatrixBuilder::addSPolynomialToMatrix
(const Poly& polyA, const Poly& polyB) {
  MATHICGB_ASSERT(!polyA.isZero());
  MATHICGB_ASSERT(polyA.isMonic());
  MATHICGB_ASSERT(!polyB.isZero());
  MATHICGB_ASSERT(polyB.isMonic());

  monomial lcm = ring().allocMonomial();
  ring().monomialLeastCommonMultiple
    (polyA.getLeadMonomial(), polyB.getLeadMonomial(), lcm);

  SPairTask task = {};

  task.polyA = &polyA;
  task.multiplyA = ring().allocMonomial();
  ring().monomialDivide(lcm, polyA.getLeadMonomial(), task.multiplyA);

  task.polyB = &polyB;
  task.multiplyB = ring().allocMonomial();
  ring().monomialDivide(lcm, polyB.getLeadMonomial(), task.multiplyB);

  mSPairTodo.push_back(task);
  ring().freeMonomial(lcm);
}

void F4MatrixBuilder::addPolynomialToMatrix(const Poly& poly) {
  if (poly.isZero())
    return;

  SPairTask task = {};

  task.polyA = &poly;
  task.multiplyA = ring().allocMonomial();
  ring().monomialSetIdentity(task.multiplyA);

  MATHICGB_ASSERT(task.polyB == 0);
  MATHICGB_ASSERT(task.multiplyB.unsafeGetRepresentation() == 0);
  mSPairTodo.push_back(task);
}

void F4MatrixBuilder::addPolynomialToMatrix
(const_monomial multiple, const Poly& poly) {
  MATHICGB_ASSERT(ring().hashValid(multiple));
  if (poly.isZero())
    return;

  SPairTask task = {};

  task.polyA = &poly;
  task.multiplyA = ring().allocMonomial();
  ring().monomialCopy(multiple, task.multiplyA);

  MATHICGB_ASSERT(task.polyB == 0);
  MATHICGB_ASSERT(task.multiplyB.unsafeGetRepresentation() == 0);
  mSPairTodo.push_back(task);
}

void F4MatrixBuilder::buildMatrixAndClear(QuadMatrix& matrix) {
  // todo: detect and remove duplicate input rows.
  // todo: prefer sparse/old reducers among the inputs.
  monomial mono = ring().allocMonomial();

  for (auto it = mSPairTodo.begin(); it != mSPairTodo.end(); ++it) {
    Poly::const_iterator itA = it->polyA->begin();
    Poly::const_iterator endA = it->polyA->end();
    MATHICGB_ASSERT(itA != endA);
    MATHICGB_ASSERT(it->multiplyA.unsafeGetRepresentation() != 0);
    MATHICGB_ASSERT(ring().hashValid(it->multiplyA));

    // make itB==endB if there is no polyB
    Poly::const_iterator itB;
    Poly::const_iterator endB;
    if (it->polyB != 0) {
      MATHICGB_ASSERT(it->multiplyB.unsafeGetRepresentation() != 0);
      MATHICGB_ASSERT(ring().hashValid(it->multiplyB));
      itB = it->polyB->begin();
      endB = it->polyB->end();
      MATHICGB_ASSERT(itB != endB);

      // skip leading terms since they cancel
      ++itA;
      ++itB;
    } else {
      // set polyB to an empty range if null
      itB = endA;
      endB = endA;
      MATHICGB_ASSERT(itB == endB);
    }

    monomial monoA = ring().allocMonomial();
    monomial monoB = ring().allocMonomial();
    monomial mono = ring().allocMonomial();

    if (itA != endA)
      ring().monomialMult(itA.getMonomial(), it->multiplyA, monoA);
    if (itB != endB)
      ring().monomialMult(itB.getMonomial(), it->multiplyB, monoB);

    while (true) {
      bool popA = false;
      bool popB = false;
      if (itA == endA) {
        if (itB == endB)
          break;
        popB = true;
      } else if (itB == endB)
        popA = true;
      else {
        int cmp = ring().monomialCompare(monoA, monoB);
        if (cmp != GT)
          popB = true;
        if (cmp != LT)
          popA = true;
      }

      MATHICGB_ASSERT(popA || popB);
      coefficient c = 0;
      if (popA) {
        std::swap(mono, monoA);
        ring().coefficientAddTo(c, itA.getCoefficient());
        ++itA;
        if (itA != endA)
          ring().monomialMult(itA.getMonomial(), it->multiplyA, monoA);
      }
      if (popB) {
        std::swap(mono, monoB);
        coefficient cB = itB.getCoefficient();
        ring().coefficientNegateTo(cB);
        ring().coefficientAddTo(c, cB);
        ++itB;
        if (itB != endB)
          ring().monomialMult(itB.getMonomial(), it->multiplyB, monoB);
      }
      if (c != 0) {
	    MATHICGB_ASSERT(c < std::numeric_limits<QuadMatrixBuilder::Scalar>::max());
        mBuilder.appendEntryBottom(createOrFindColumnOf(mono), static_cast<QuadMatrixBuilder::Scalar>(c));
	  }
    }
    ring().freeMonomial(monoA);
    ring().freeMonomial(monoB);
    mBuilder.rowDoneBottomLeftAndRight();
  }

  // Process pending rows until we are done. Note that the methods
  // we are calling here can add more items to mTodo.
  while (!mTodo.empty()) {
    RowTask task = mTodo.back();
    MATHICGB_ASSERT(ring().hashValid(task.multiple));
    mTodo.pop_back();
    if (task.useAsReducer)
      appendRowTop(task.multiple, *task.poly);
    else
      appendRowBottom(task.multiple, *task.poly);
  }

  mBuilder.sortColumnsLeft(mBasis.order());
  mBuilder.sortColumnsRight(mBasis.order());
  mBuilder.buildMatrixAndClear(matrix);
}

F4MatrixBuilder::LeftRightColIndex
F4MatrixBuilder::createOrFindColumnOf(const_monomial mono) {
  MATHICGB_ASSERT(ring().hashValid(mono));
  LeftRightColIndex colIndex = mBuilder.findColumn(mono);
  if (colIndex.valid())
    return colIndex;

  // mono did not already have a column so look for a reducer
  size_t reducerIndex = mBasis.divisor(mono);
  if (reducerIndex == static_cast<size_t>(-1))
    return mBuilder.createColumnRight(mono);

  // schedule the reducer to be added as a row
  RowTask task;
  task.poly = &mBasis.poly(reducerIndex);
  task.useAsReducer = true;
  task.multiple = ring().allocMonomial();
  ring().monomialDivideToNegative
    (mono, task.poly->getLeadMonomial(), task.multiple);
  MATHICGB_ASSERT(ring().hashValid(task.multiple));
  mTodo.push_back(task);

  return mBuilder.createColumnLeft(mono);
}

void F4MatrixBuilder::appendRowTop(const_monomial multiple, const Poly& poly) {
  monomial mono = ring().allocMonomial();
  Poly::const_iterator end = poly.end();
  for (Poly::const_iterator it = poly.begin(); it != end; ++it) {
    ring().monomialMult(it.getMonomial(), multiple, mono);
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<QuadMatrixBuilder::Scalar>::max());
    mBuilder.appendEntryTop(createOrFindColumnOf(mono), static_cast<QuadMatrixBuilder::Scalar>(it.getCoefficient()));
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
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<QuadMatrixBuilder::Scalar>::max());
    mBuilder.appendEntryBottom
      (createOrFindColumnOf(mono), static_cast<QuadMatrixBuilder::Scalar>(it.getCoefficient()));
  }
  mBuilder.rowDoneBottomLeftAndRight();
  ring().freeMonomial(mono);
}
