#include "stdinc.h"
#include "F4MatrixBuilder.hpp"

#include <tbb/tbb.h>

MATHICGB_NO_INLINE
std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial>
F4MatrixBuilder::findOrCreateColumn(
  const const_monomial monoA,
  const const_monomial monoB
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());
  const auto col = ColReader(mMap).findProduct(monoA, monoB);
  if (col.first != 0)
    return std::make_pair(*col.first, col.second);
  return createColumn(monoA, monoB);
}

MATHICGB_INLINE
std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial>
F4MatrixBuilder::findOrCreateColumn(
  const const_monomial monoA,
  const const_monomial monoB,
  const ColReader& colMap
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());
  const auto col = colMap.findProduct(monoA, monoB);
  if (col.first == 0)
    return findOrCreateColumn(monoA, monoB);
  return std::make_pair(*col.first, col.second);
}

MATHICGB_NO_INLINE
void F4MatrixBuilder::createTwoColumns(
  const const_monomial monoA1,
  const const_monomial monoA2,
  const const_monomial monoB
) {
  createColumn(monoA1, monoB);
  createColumn(monoA2, monoB);
}

F4MatrixBuilder::F4MatrixBuilder(
  const PolyBasis& basis,
  const int threadCount,
  const size_t memoryQuantum
):
  mThreadCount(threadCount),
  mTmp(basis.ring().allocMonomial()),
  mBasis(basis),
  mMap(basis.ring()),
  mMonomialsLeft(),
  mMonomialsRight(),
  mBuilder(basis.ring(), mMap, mMonomialsLeft, mMonomialsRight, memoryQuantum),
  mLeftColCount(0),
  mRightColCount(0)
{
  MATHICGB_ASSERT(threadCount >= 1);

  // This assert to be _NO_ASSUME since otherwise the compiler will assume that
  // the error checking branch here cannot be taken and optimize it away.
  const Scalar maxScalar = std::numeric_limits<Scalar>::max();
  MATHICGB_ASSERT_NO_ASSUME(ring().charac() <= maxScalar);
  if (ring().charac() > maxScalar)
    mathic::reportInternalError("F4MatrixBuilder: too large characteristic.");
}

void F4MatrixBuilder::addSPolynomialToMatrix
(const Poly& polyA, const Poly& polyB) {
  MATHICGB_ASSERT(!polyA.isZero());
  MATHICGB_ASSERT(polyA.isMonic());
  MATHICGB_ASSERT(!polyB.isZero());
  MATHICGB_ASSERT(polyB.isMonic());

  monomial lcm = ring().allocMonomial();
  ring().monomialLeastCommonMultiple
    (polyA.getLeadMonomial(), polyB.getLeadMonomial(), lcm);

  RowTask task;
  task.addToTop = false;

  task.poly = &polyA;
  task.multiply = ring().allocMonomial();
  ring().monomialDivide(lcm, polyA.getLeadMonomial(), task.multiply);

  task.sPairPoly = &polyB;
  task.sPairMultiply = ring().allocMonomial();
  ring().monomialDivide(lcm, polyB.getLeadMonomial(), task.sPairMultiply);

  mTodo.push_back(task);
  ring().freeMonomial(lcm);
}

void F4MatrixBuilder::addPolynomialToMatrix(const Poly& poly) {
  if (poly.isZero())
    return;

  RowTask task = {};
  task.addToTop = false;

  task.poly = &poly;
  task.multiply = ring().allocMonomial();
  ring().monomialSetIdentity(task.multiply);

  MATHICGB_ASSERT(task.sPairPoly == 0);
  MATHICGB_ASSERT(task.sPairMultiply.isNull());
  mTodo.push_back(task);
}

void F4MatrixBuilder::addPolynomialToMatrix
(const_monomial multiple, const Poly& poly) {
  MATHICGB_ASSERT(ring().hashValid(multiple));
  if (poly.isZero())
    return;

  RowTask task = {};
  task.addToTop = false;

  task.poly = &poly;
  task.multiply = ring().allocMonomial();
  ring().monomialCopy(multiple, task.multiply);

  MATHICGB_ASSERT(task.sPairPoly == 0);
  MATHICGB_ASSERT(task.sPairMultiply.isNull());
  mTodo.push_back(task);
}

void F4MatrixBuilder::buildMatrixAndClear(QuadMatrix& matrix) {
  // todo: prefer sparse/old reducers among the inputs.

  // Process pending rows until we are done. Note that the methods
  // we are calling here can add more pending items.

  MATHICGB_ASSERT(mThreadCount >= 1);
  tbb::enumerable_thread_specific<std::unique_ptr<QuadMatrixBuilder>>
  builders([&](){
    return make_unique<QuadMatrixBuilder>(
      ring(), mMap, mMonomialsLeft, mMonomialsRight, mBuilder.memoryQuantum()
    );
  }); 

  decltype(mTodo) currentTasks;
  while (!mTodo.empty()) {
    for (auto it = currentTasks.begin(); it != currentTasks.end(); ++it) {
      MATHICGB_ASSERT(!it->multiply.isNull());
      ring().freeMonomial(it->multiply);
      MATHICGB_ASSERT(it->sPairMultiply.isNull() == (it->sPairPoly == 0));
      if (!it->sPairMultiply.isNull())
        ring().freeMonomial(it->sPairMultiply);
    }
    currentTasks.clear();
    mTodo.swap(currentTasks);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, currentTasks.size(), 1),
      [&](const tbb::blocked_range<size_t>& range)
      {for (auto it = range.begin(); it != range.end(); ++it)
    {
      const size_t taskIndex = it;
      QuadMatrixBuilder& builder = *builders.local();
      MATHICGB_ASSERT(&builder != 0);

      const RowTask task = currentTasks[taskIndex];
      MATHICGB_ASSERT(ring().hashValid(task.multiply));

      if (task.addToTop) {
        MATHICGB_ASSERT(task.sPairPoly == 0);
        appendRowTop(task.multiply, *task.poly, builder);
      } else
        appendRowBottom(task, builder);
    }});
  }

  if (builders.empty()) {
    matrix = QuadMatrix();
    matrix.ring = &ring();
    return;
  }

  auto& builder = **builders.begin();
  const auto end = builders.end();
  for (auto it = builders.begin() + 1; it != end; ++it)
    builder.takeRowsFrom((*it)->buildMatrixAndClear());
  matrix = builder.buildMatrixAndClear();
  builders.clear();

  {
    ColReader reader(mMap);
    matrix.leftColumnMonomials.clear();
    matrix.rightColumnMonomials.clear();
    const auto end = reader.end();
    for (auto it = reader.begin(); it != end; ++it) {
      const auto p = *it;
      monomial copy = ring().allocMonomial();
      ring().monomialCopy(p.second, copy);
      auto& monos = p.first.left() ?
        matrix.leftColumnMonomials : matrix.rightColumnMonomials;
      const auto index = p.first.index();
      if (monos.size() <= index)
        monos.resize(index + 1);
      MATHICGB_ASSERT(monos[index].isNull());
      monos[index] = copy;
    }
  }
#ifdef MATHICGB_DEBUG
  for (size_t side = 0; side < 2; ++side) {
    auto& monos = side == 0 ?
      matrix.leftColumnMonomials : matrix.rightColumnMonomials;
    for (auto it = monos.begin(); it != monos.end(); ++it) {
      MATHICGB_ASSERT(!it->isNull());
    }
  }
#endif
  matrix.sortColumnsLeftRightParallel(mThreadCount);
  mMap.clearNonConcurrent();
}

std::pair<F4MatrixBuilder::LeftRightColIndex, ConstMonomial>
F4MatrixBuilder::createColumn(
  const const_monomial monoA,
  const const_monomial monoB
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());

  tbb::mutex::scoped_lock lock(mCreateColumnLock);
  // see if the column exists now after we have synchronized
  {
    const auto found(ColReader(mMap).findProduct(monoA, monoB));
    if (found.first != 0)
      return std::make_pair(*found.first, found.second);
  }

  // The column really does not exist, so we need to create it
  ring().monomialMult(monoA, monoB, mTmp);
  MATHICGB_ASSERT(ring().hashValid(mTmp));

  // look for a reducer of mTmp
  const size_t reducerIndex = mBasis.divisor(mTmp);
  const bool insertLeft = (reducerIndex != static_cast<size_t>(-1));
  if (insertLeft) {
    RowTask task = {}; // schedule the new reducer to be added as a row
    task.addToTop = true;
    task.poly = &mBasis.poly(reducerIndex);
    task.multiply = ring().allocMonomial();
    ring().monomialDivideToNegative
      (mTmp, task.poly->getLeadMonomial(), task.multiply);
    MATHICGB_ASSERT(ring().hashValid(task.multiply));
    mTodo.push_back(task);
  }

  // Create the new left or right column
  const auto colCount = insertLeft ? mLeftColCount : mRightColCount;
  if (colCount == std::numeric_limits<ColIndex>::max())
    throw std::overflow_error("Too many columns in QuadMatrix");
  const auto inserted = mMap.insert
    (std::make_pair(mTmp, LeftRightColIndex(colCount, insertLeft)));
  insertLeft ? ++mLeftColCount : ++mRightColCount;
  MATHICGB_ASSERT(inserted.second);
  MATHICGB_ASSERT(inserted.first.first != 0);
  return std::make_pair(*inserted.first.first, inserted.first.second);
}

void F4MatrixBuilder::appendRowBottom(
  const_monomial multiple,
  const bool negate,
  const Poly::const_iterator begin,
  const Poly::const_iterator end,
  QuadMatrixBuilder& builder
) {
  // todo: eliminate the code-duplication between here and appendRowTop.
  MATHICGB_ASSERT(!multiple.isNull());
  MATHICGB_ASSERT(&builder != 0);

  auto it = begin;
updateReader:
  // Use an on-stack const reader to make it as obvious as possible to the
  // optimizer's alias analysis that the pointer inside the reader never
  // changes inside the loop.
  const ColReader reader(mMap);
  for (; it != end; ++it) {
    const auto col = reader.findProduct(it.getMonomial(), multiple);
    if (col.first == 0) {
      createColumn(it.getMonomial(), multiple);
      goto updateReader;
    }

    const auto origScalar = it.getCoefficient();
    MATHICGB_ASSERT(origScalar != 0);
    const auto maybeNegated =
      negate ? ring().coefficientNegateNonZero(origScalar) : origScalar;
	MATHICGB_ASSERT(maybeNegated < std::numeric_limits<Scalar>::max());
    builder.appendEntryBottom(*col.first, static_cast<Scalar>(maybeNegated));
  }
  builder.rowDoneBottomLeftAndRight();
}

void F4MatrixBuilder::appendRowTop(
  const const_monomial multiple,
  const Poly& poly,
  QuadMatrixBuilder& builder
) {
  MATHICGB_ASSERT(!multiple.isNull());
  MATHICGB_ASSERT(&poly != 0);
  MATHICGB_ASSERT(&builder != 0);

  auto it = poly.begin();
  const auto end = poly.end();
  if ((std::distance(it, end) % 2) == 1) {
    ColReader reader(mMap);
    const auto col = findOrCreateColumn(it.getMonomial(), multiple, reader);
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient());
    builder.appendEntryTop(col.first, static_cast<Scalar>(it.getCoefficient()));
    ++it;
  }
updateReader:
  ColReader colMap(mMap);
  MATHICGB_ASSERT((std::distance(it, end) % 2) == 0);
  while (it != end) {
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient() != 0);
    const auto scalar1 = static_cast<Scalar>(it.getCoefficient());
    const const_monomial mono1 = it.getMonomial();

    auto it2 = it;
    ++it2;
	MATHICGB_ASSERT(it2.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it2.getCoefficient() != 0);
    const auto scalar2 = static_cast<Scalar>(it2.getCoefficient());
    const const_monomial mono2 = it2.getMonomial();

    const auto colPair = colMap.findTwoProducts(mono1, mono2, multiple);
    if (colPair.first == 0 || colPair.second == 0) {
      createTwoColumns(mono1, mono2, multiple);
      goto updateReader;
    }

    builder.appendEntryTop(*colPair.first, scalar1);
    builder.appendEntryTop(*colPair.second, scalar2);
    it = ++it2;
  }
  builder.rowDoneTopLeftAndRight();
}

void F4MatrixBuilder::appendRowBottom(
  const RowTask& task,
  QuadMatrixBuilder& builder
) {
  MATHICGB_ASSERT(!task.addToTop);
  MATHICGB_ASSERT(!task.poly->isZero());
  MATHICGB_ASSERT(!task.multiply.isNull());
  MATHICGB_ASSERT(ring().hashValid(task.multiply));
  Poly::const_iterator itA = task.poly->begin();
  const Poly::const_iterator endA = task.poly->end();

  if (task.sPairPoly == 0) {
    appendRowBottom(task.multiply, false, itA, endA, builder);
    return;
  }

  MATHICGB_ASSERT(!task.sPairPoly->isZero());
  MATHICGB_ASSERT(!task.sPairMultiply.isNull());
  MATHICGB_ASSERT(ring().hashValid(task.sPairMultiply));
  Poly::const_iterator itB = task.sPairPoly->begin();
  Poly::const_iterator endB = task.sPairPoly->end();

  // skip leading terms since they cancel
  MATHICGB_ASSERT(itA.getCoefficient() == itB.getCoefficient());
  ++itA;
  ++itB;

  const ColReader colMap(mMap);

  const const_monomial mulA = task.multiply;
  const const_monomial mulB = task.sPairMultiply;
  while (true) {
    // Watch out: we are depending on appendRowBottom to finish the row, so
    // if you decide not to call that function in case
    // (itA == itA && itB == endB) then you need to do that yourself.
    if (itB == endB) {
      appendRowBottom(mulA, false, itA, endA, builder);
      break;
    }
    if (itA == endA) {
      appendRowBottom(mulB, true, itB, endB, builder);
      break;
    }

    coefficient coeff = 0;
    LeftRightColIndex col;
    const auto colA = findOrCreateColumn(itA.getMonomial(), mulA, colMap);
    const auto colB = findOrCreateColumn(itB.getMonomial(), mulB, colMap);
    const auto cmp = ring().monomialCompare(colA.second, colB.second);
    //const auto cmp = ring().monomialCompare
    //  (builder.monomialOfCol(colA), builder.monomialOfCol(colB));
    if (cmp != LT) {
      coeff = itA.getCoefficient();
      col = colA.first;
      ++itA;
    }
    if (cmp != GT) {
      coeff = ring().coefficientSubtract(coeff, itB.getCoefficient());
      col = colB.first;
      ++itB;
    }
    MATHICGB_ASSERT(coeff < std::numeric_limits<Scalar>::max());
    if (coeff != 0)
      builder.appendEntryBottom(col, static_cast<Scalar>(coeff));
  }
}
