// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "F4MatrixBuilder.hpp"

#include "LogDomain.hpp"

MATHICGB_DEFINE_LOG_DOMAIN(
  F4MatrixBuild,
  "Displays statistics about F4 matrix construction."
);

MATHICGB_NAMESPACE_BEGIN

MATHICGB_NO_INLINE
std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial>
F4MatrixBuilder::findOrCreateColumn(
  const const_monomial monoA,
  const const_monomial monoB,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());
  const auto col = ColReader(mMap).findProduct(monoA, monoB);
  if (col.first != 0)
    return std::make_pair(*col.first, col.second);
  return createColumn(monoA, monoB, feeder);
}

MATHICGB_INLINE
std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial>
F4MatrixBuilder::findOrCreateColumn(
  const const_monomial monoA,
  const const_monomial monoB,
  const ColReader& colMap,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());
  const auto col = colMap.findProduct(monoA, monoB);
  if (col.first == 0)
    return findOrCreateColumn(monoA, monoB, feeder);
  return std::make_pair(*col.first, col.second);
}

MATHICGB_NO_INLINE
void F4MatrixBuilder::createTwoColumns(
  const const_monomial monoA1,
  const const_monomial monoA2,
  const const_monomial monoB,
  TaskFeeder& feeder
) {
  createColumn(monoA1, monoB, feeder);
  createColumn(monoA2, monoB, feeder);
}

F4MatrixBuilder::F4MatrixBuilder(
  const PolyBasis& basis,
  const size_t memoryQuantum
):
  mTmp(basis.ring().allocMonomial()),
  mBasis(basis),
  mMap(basis.ring()),
  mMonomialsLeft(),
  mMonomialsRight(),
  mBuilder(basis.ring(), mMap, mMonomialsLeft, mMonomialsRight, memoryQuantum),
  mLeftColCount(0),
  mRightColCount(0)
{
  // This assert to be _NO_ASSUME since otherwise the compiler will assume that
  // the error checking branch here cannot be taken and optimize it away.
  const Scalar maxScalar = std::numeric_limits<Scalar>::max();
  MATHICGB_ASSERT_NO_ASSUME(ring().charac() <= maxScalar);
  if (ring().charac() > maxScalar)
    mathic::reportInternalError("F4MatrixBuilder: too large characteristic.");
}

void F4MatrixBuilder::addSPolynomialToMatrix(
  const Poly& polyA,
  const Poly& polyB
) {
  MATHICGB_ASSERT(!polyA.isZero());
  MATHICGB_ASSERT(polyA.isMonic());
  MATHICGB_ASSERT(!polyB.isZero());
  MATHICGB_ASSERT(polyB.isMonic());

  RowTask task;
  task.addToTop = false;
  task.poly = &polyA;
  task.sPairPoly = &polyB;
  mTodo.push_back(task);
}

void F4MatrixBuilder::addPolynomialToMatrix(const Poly& poly) {
  if (poly.isZero())
    return;

  RowTask task = {};
  task.addToTop = false;
  task.poly = &poly;
  mTodo.push_back(task);
}

void F4MatrixBuilder::addPolynomialToMatrix
(const_monomial multiple, const Poly& poly) {
  if (poly.isZero())
    return;

  RowTask task = {};
  task.addToTop = false;
  task.poly = &poly;
  task.desiredLead = ring().allocMonomial();
  ring().monomialMult(poly.getLeadMonomial(), multiple, task.desiredLead);

  MATHICGB_ASSERT(task.sPairPoly == 0);
  mTodo.push_back(task);
}

void F4MatrixBuilder::buildMatrixAndClear(QuadMatrix& matrix) {
  MATHICGB_LOG_TIME(F4MatrixBuild) <<
    "\n***** Constructing matrix *****\n";

  if (mTodo.empty()) {
    matrix = QuadMatrix();
    matrix.ring = &ring();
    return;
  }

  // todo: prefer sparse/old reducers among the inputs.

  // Process pending rows until we are done. Note that the methods
  // we are calling here can add more pending items.

  struct ThreadData {
    QuadMatrixBuilder builder;
    monomial tmp1;
    monomial tmp2;
  };

  mgb::mtbb::enumerable_thread_specific<ThreadData> threadData([&](){  
    ThreadData data = {QuadMatrixBuilder(
      ring(), mMap, mMonomialsLeft, mMonomialsRight, mBuilder.memoryQuantum()
    )};
    {
      mgb::mtbb::mutex::scoped_lock guard(mCreateColumnLock);
      data.tmp1 = ring().allocMonomial();
      data.tmp2 = ring().allocMonomial();
    }
    return std::move(data);
  });

  mgb::mtbb::parallel_do(mTodo.begin(), mTodo.end(),
    [&](const RowTask& task, TaskFeeder& feeder)
  {
    auto& data = threadData.local();
    QuadMatrixBuilder& builder = data.builder;
    const Poly& poly = *task.poly;

    if (task.sPairPoly != 0) {
      MATHICGB_ASSERT(!task.addToTop);
      ring().monomialColons(
        poly.getLeadMonomial(),
        task.sPairPoly->getLeadMonomial(),
        data.tmp2,
        data.tmp1
      );
      appendRowBottom
        (&poly, data.tmp1, task.sPairPoly, data.tmp2, data.builder, feeder);
      return;
    }
    if (task.desiredLead.isNull())
      ring().monomialSetIdentity(data.tmp1);
    else
      ring().monomialDivide
        (task.desiredLead, poly.getLeadMonomial(), data.tmp1);
    if (task.addToTop)
      appendRowTop(data.tmp1, *task.poly, builder, feeder);
    else
      appendRowBottom
        (data.tmp1, false, poly.begin(), poly.end(), builder, feeder);
  });
  MATHICGB_ASSERT(!threadData.empty()); // as mTodo empty causes early return
  // Free the monomials from all the tasks
  const auto todoEnd = mTodo.end();
  for (auto it = mTodo.begin(); it != todoEnd; ++it)
    if (!it->desiredLead.isNull())
      ring().freeMonomial(it->desiredLead);
  mTodo.clear();

  auto& builder = threadData.begin()->builder;
  const auto end = threadData.end();
  for (auto it = threadData.begin(); it != end; ++it) {
    if (&it->builder != &builder)
      builder.takeRowsFrom(it->builder.buildMatrixAndClear());
    ring().freeMonomial(it->tmp1);
    ring().freeMonomial(it->tmp2);
  }
  matrix = builder.buildMatrixAndClear();
  threadData.clear();

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
  matrix.sortColumnsLeftRightParallel();
  mMap.clearNonConcurrent();
}

std::pair<F4MatrixBuilder::LeftRightColIndex, ConstMonomial>
F4MatrixBuilder::createColumn(
  const const_monomial monoA,
  const const_monomial monoB,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());

  mgb::mtbb::mutex::scoped_lock lock(mCreateColumnLock);
  // see if the column exists now after we have synchronized
  {
    const auto found(ColReader(mMap).findProduct(monoA, monoB));
    if (found.first != 0)
      return std::make_pair(*found.first, found.second);
  }

  // The column really does not exist, so we need to create it
  ring().monomialMult(monoA, monoB, mTmp);
  if (!ring().monomialHasAmpleCapacity(mTmp))
    mathic::reportError("Monomial exponent overflow in F4MatrixBuilder.");

  // look for a reducer of mTmp
  const size_t reducerIndex = mBasis.classicReducer(mTmp);
  const bool insertLeft = (reducerIndex != static_cast<size_t>(-1));

  // Create the new left or right column
  auto& colCount = insertLeft ? mLeftColCount : mRightColCount;
  if (colCount == std::numeric_limits<ColIndex>::max())
    throw std::overflow_error("Too many columns in QuadMatrix");
  const auto inserted = mMap.insert
    (std::make_pair(mTmp, LeftRightColIndex(colCount, insertLeft)));
  ++colCount;
  MATHICGB_ASSERT(inserted.second);
  MATHICGB_ASSERT(inserted.first.first != 0);

  // schedule new task if we found a reducer
  if (insertLeft) {
    RowTask task = {};
    task.addToTop = true;
    task.poly = &mBasis.poly(reducerIndex);
    task.desiredLead = inserted.first.second.castAwayConst();
    feeder.add(task);
  }

  return std::make_pair(*inserted.first.first, inserted.first.second);
}

void F4MatrixBuilder::appendRowBottom(
  const_monomial multiple,
  const bool negate,
  const Poly::const_iterator begin,
  const Poly::const_iterator end,
  QuadMatrixBuilder& builder,
  TaskFeeder& feeder
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
      createColumn(it.getMonomial(), multiple, feeder);
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
  QuadMatrixBuilder& builder,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!multiple.isNull());
  MATHICGB_ASSERT(&poly != 0);
  MATHICGB_ASSERT(&builder != 0);

  auto it = poly.begin();
  const auto end = poly.end();
  if ((std::distance(it, end) % 2) == 1) {
    ColReader reader(mMap);
    const auto col = findOrCreateColumn
      (it.getMonomial(), multiple, reader, feeder);
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient());
    builder.appendEntryTop
      (col.first, static_cast<Scalar>(it.getCoefficient()));
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
      createTwoColumns(mono1, mono2, multiple, feeder);
      goto updateReader;
    }

    builder.appendEntryTop(*colPair.first, scalar1);
    builder.appendEntryTop(*colPair.second, scalar2);
    it = ++it2;
  }
  builder.rowDoneTopLeftAndRight();
}

void F4MatrixBuilder::appendRowBottom(
  const Poly* poly,
  monomial multiply,
  const Poly* sPairPoly,
  monomial sPairMultiply,
  QuadMatrixBuilder& builder,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!poly->isZero());
  MATHICGB_ASSERT(!multiply.isNull());
  MATHICGB_ASSERT(sPairPoly != 0);
  Poly::const_iterator itA = poly->begin();
  const Poly::const_iterator endA = poly->end();

  MATHICGB_ASSERT(!sPairPoly->isZero());
  MATHICGB_ASSERT(!sPairMultiply.isNull());
  Poly::const_iterator itB = sPairPoly->begin();
  Poly::const_iterator endB = sPairPoly->end();

  // skip leading terms since they cancel
  MATHICGB_ASSERT(itA.getCoefficient() == itB.getCoefficient());
  ++itA;
  ++itB;

  const ColReader colMap(mMap);

  const const_monomial mulA = multiply;
  const const_monomial mulB = sPairMultiply;
  while (true) {
    // Watch out: we are depending on appendRowBottom to finish the row, so
    // if you decide not to call that function in case
    // (itA == itA && itB == endB) then you need to do that yourself.
    if (itB == endB) {
      appendRowBottom(mulA, false, itA, endA, builder, feeder);
      break;
    }
    if (itA == endA) {
      appendRowBottom(mulB, true, itB, endB, builder, feeder);
      break;
    }

    coefficient coeff = 0;
    LeftRightColIndex col;
    const auto colA = findOrCreateColumn
      (itA.getMonomial(), mulA, colMap, feeder);
    const auto colB = findOrCreateColumn
      (itB.getMonomial(), mulB, colMap, feeder);
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

MATHICGB_NAMESPACE_END
