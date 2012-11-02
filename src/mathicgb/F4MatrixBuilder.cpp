#include "stdinc.h"
#include "F4MatrixBuilder.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

MATHICGB_INLINE QuadMatrixBuilder::LeftRightColIndex
  F4MatrixBuilder::findOrCreateColumn
(
  const const_monomial monoA,
  const const_monomial monoB,
  QuadMatrixBuilder& builder
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());
  const auto col = builder.findColumnProduct(monoA, monoB);
  if (col.valid())
    return col;
  return createColumn(builder, monoA, monoB);
}

MATHICGB_INLINE const std::pair<
  QuadMatrixBuilder::LeftRightColIndex,
  QuadMatrixBuilder::LeftRightColIndex
> F4MatrixBuilder::findOrCreateTwoColumns
(
  const const_monomial monoA1,
  const const_monomial monoA2,
  const const_monomial monoB,
  QuadMatrixBuilder& builder
) {
  MATHICGB_ASSERT(!monoA1.isNull());
  MATHICGB_ASSERT(!monoA2.isNull());
  MATHICGB_ASSERT(!monoB.isNull());
  MATHICGB_ASSERT(!ring().monomialEQ(monoA1, monoA2));
  auto colPair = builder.findTwoColumnsProduct(monoA1, monoA2, monoB);
  if (!colPair.first.valid()) {
    colPair.first = createColumn(builder, monoA1, monoB);
    // Syncing builder to mBuilder could have created a col for monoA2*monoB.
    if (&mBuilder != &builder && !colPair.second.valid())
      colPair.second = builder.findColumnProduct(monoA2, monoB);
  }
  if (!colPair.second.valid())
    colPair.second = createColumn(builder, monoA2, monoB);
  MATHICGB_ASSERT(colPair == std::make_pair(
    findOrCreateColumn(monoA1, monoB, builder),
    findOrCreateColumn(monoA2, monoB, builder)));
  return colPair;
}

F4MatrixBuilder::F4MatrixBuilder(
  const PolyBasis& basis,
  const int threadCount,
  const size_t memoryQuantum
):
  mThreadCount(threadCount),
  mBasis(basis),
  mBuilder(basis.ring(), memoryQuantum),
  mTmp(basis.ring().allocMonomial())
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
  // we are calling here can add more items to mTodo.

#ifdef _OPENMP
  struct ThreadData {
    QuadMatrixBuilder* builder;
    monomial tmp;
  };
  MATHICGB_ASSERT(mThreadCount >= 1);
  std::vector<ThreadData> threadData(mThreadCount);
  if (mThreadCount == 1) {
    threadData[0].builder = &mBuilder;
    threadData[0].tmp = mTmp;
  } else {
    for (size_t i = 0; i < threadData.size(); ++i) {
      threadData[i].builder = new QuadMatrixBuilder(ring());
      threadData[i].tmp = new exponent[ring().maxMonomialByteSize()];
    }
  }
#endif

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
    const auto taskCountOMP = static_cast<OMPIndex>(currentTasks.size());
#pragma omp parallel for num_threads(mThreadCount) schedule(dynamic)
    for (OMPIndex taskOMP = 0; taskOMP < taskCountOMP; ++taskOMP) {
#pragma omp flush

      const size_t taskIndex = taskOMP;
#ifdef _OPENMP
      MATHICGB_ASSERT(omp_get_thread_num() < threadData.size());
      const ThreadData& td = threadData[omp_get_thread_num()];
      QuadMatrixBuilder& builder = *td.builder;
#else
      QuadMatrixBuilder& builder = mBuilder;
#endif

      const RowTask task = currentTasks[taskIndex];
      MATHICGB_ASSERT(ring().hashValid(task.multiply));

      if (task.addToTop) {
        MATHICGB_ASSERT(task.sPairPoly == 0);
        appendRowTop(task.multiply, *task.poly, builder);
      } else
        appendRowBottom(task, builder);
    }
  }

#ifdef _OPENMP
#pragma omp flush
  MATHICGB_ASSERT(mThreadCount > 1 || threadData.back().builder == &mBuilder);
  if (mThreadCount > 1) {
    for (auto it = threadData.begin(); it != threadData.end(); ++it) {
      MATHICGB_ASSERT(&mBuilder != it->builder);
      QuadMatrix qm(it->builder->buildMatrixAndClear());
      mBuilder.takeRowsFrom(std::move(qm));
      delete it->builder;
    }
  }
  threadData.clear();
#endif

  mBuilder.sortColumnsLeftRightParallel(mBasis.order(), mThreadCount);
  matrix = mBuilder.buildMatrixAndClear();
}

F4MatrixBuilder::LeftRightColIndex
F4MatrixBuilder::createColumn(
  QuadMatrixBuilder& builder,
  const const_monomial monoA,
  const const_monomial monoB
) {
  MATHICGB_ASSERT(!monoA.isNull());
  MATHICGB_ASSERT(!monoB.isNull());
  // Must be declared up here since we cannot return out of an omp critical
  // section.
  QuadMatrixBuilder::LeftRightColIndex newCol;

#pragma omp critical
  {
    ring().monomialMult(monoA, monoB, mTmp);

    MATHICGB_ASSERT(ring().hashValid(mTmp));
    MATHICGB_ASSERT(!builder.findColumn(mTmp).valid());
#ifndef _OPENMP
    MATHICGB_ASSERT(&builder == &mBuilder);
#endif

#ifdef _OPENMP
    // builder does not have a column for mTmp but mBuilder might
    if (&mBuilder != &builder)
      newCol = mBuilder.findColumn(mTmp);
    if (!newCol.valid())
#endif
    {
      // we need to add a new column to mBuilder

      // look for a reducer of mTmp
      size_t reducerIndex = mBasis.divisor(mTmp);
      if (reducerIndex == static_cast<size_t>(-1)) {
        newCol = mBuilder.createColumnRight(mTmp);
      } else {
        newCol = mBuilder.createColumnLeft(mTmp);

        // schedule the reducer to be added as a row
        RowTask task = {};
        task.addToTop = true;
        task.poly = &mBasis.poly(reducerIndex);
        task.multiply = ring().allocMonomial();
        ring().monomialDivideToNegative
          (mTmp, task.poly->getLeadMonomial(), task.multiply);
        MATHICGB_ASSERT(ring().hashValid(task.multiply));
        mTodo.push_back(task);
      }
    }
    MATHICGB_ASSERT(newCol.valid());
    MATHICGB_ASSERT(ring().monomialEQ(mBuilder.monomialOfCol(newCol), mTmp));
    MATHICGB_ASSERT(newCol == mBuilder.findColumn(mTmp));

#ifdef _OPENMP
    // If we are running in parallel we need to catch builder up to date with
    // with respect to the columns in mBuilder. We have just ensured that
    // mBuilder has a column for mTmp so now builder will also get that column.

    if (&builder != &mBuilder) {
      // add missing left columns
      if (builder.leftColCount() != mBuilder.leftColCount()) {
        const auto colCount = mBuilder.leftColCount();
        MATHICGB_ASSERT(builder.leftColCount() < colCount);
        for (auto col = builder.leftColCount(); col < colCount; ++col) {
          const const_monomial monoToCopy = mBuilder.monomialOfLeftCol(col);
          MATHICGB_ASSERT(!builder.findColumn(monoToCopy).valid());
          builder.createColumnLeft(monoToCopy);
          MATHICGB_ASSERT(builder.findColumn(monoToCopy) ==
            mBuilder.findColumn(monoToCopy));
          MATHICGB_ASSERT(ring().monomialEQ
            (mBuilder.monomialOfLeftCol(col), builder.monomialOfLeftCol(col)));
        }
      }
#ifdef MATHICGB_DEBUG
      MATHICGB_ASSERT(builder.leftColCount() == mBuilder.leftColCount());
      for (SparseMatrix::ColIndex col = 0;
        col < builder.leftColCount(); ++col) {
        MATHICGB_ASSERT(ring().monomialEQ
          (mBuilder.monomialOfLeftCol(col), builder.monomialOfLeftCol(col)));
      }
#endif

      // add missing right columns
      if (builder.rightColCount() != mBuilder.rightColCount()) {
        const auto colCount = mBuilder.rightColCount();
        MATHICGB_ASSERT(builder.rightColCount() < colCount);
        for (auto col = builder.rightColCount(); col < colCount; ++col) {
          const const_monomial monoToCopy = mBuilder.monomialOfRightCol(col);
          MATHICGB_ASSERT(!builder.findColumn(monoToCopy).valid());
          builder.createColumnRight(monoToCopy);
          MATHICGB_ASSERT(builder.findColumn(monoToCopy) ==
            mBuilder.findColumn(monoToCopy));
          MATHICGB_ASSERT(ring().monomialEQ
            (mBuilder.monomialOfRightCol(col), builder.monomialOfRightCol(col)));
        }
      }
#ifdef MATHICGB_DEBUG
      MATHICGB_ASSERT(builder.rightColCount() == mBuilder.rightColCount());
      for (SparseMatrix::ColIndex col = 0;
        col < builder.rightColCount(); ++col) {
        MATHICGB_ASSERT(ring().monomialEQ
          (mBuilder.monomialOfRightCol(col), builder.monomialOfRightCol(col)));
      }
#endif

      MATHICGB_ASSERT(ring().monomialEQ(mBuilder.monomialOfCol(mBuilder.findColumn(mTmp)), mTmp));
      MATHICGB_ASSERT(newCol == mBuilder.findColumn(mTmp));

      MATHICGB_ASSERT(builder.findColumn(mTmp).valid());
      MATHICGB_ASSERT(ring().monomialEQ(builder.monomialOfCol(builder.findColumn(mTmp)), mTmp));
      MATHICGB_ASSERT(newCol == builder.findColumn(mTmp));
    }
#endif
  }
  return newCol;
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

  for (auto it  = begin; it != end; ++it) {
    const auto col = findOrCreateColumn(it.getMonomial(), multiple, builder);
    const auto origScalar = it.getCoefficient();
    MATHICGB_ASSERT(origScalar != 0);
    const auto possiblyNegated =
      negate ? ring().coefficientNegateNonZero(origScalar) : origScalar;
	MATHICGB_ASSERT(possiblyNegated < std::numeric_limits<Scalar>::max());
    builder.appendEntryBottom(col, static_cast<Scalar>(possiblyNegated));
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
    const auto col = findOrCreateColumn(it.getMonomial(), multiple, builder);
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient());
    builder.appendEntryTop(col, static_cast<Scalar>(it.getCoefficient()));
    ++it;
  }
  MATHICGB_ASSERT((std::distance(it, end) % 2) == 0);
  while (it != end) {
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient() != 0);
    const auto scalar1 = static_cast<Scalar>(it.getCoefficient());
    const const_monomial mono1 = it.getMonomial();
    ++it;

	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient() != 0);
    const auto scalar2 = static_cast<Scalar>(it.getCoefficient());
    const const_monomial mono2 = it.getMonomial();
    ++it;

    const auto colPair =
      findOrCreateTwoColumns(mono1, mono2, multiple, builder);
    builder.appendEntryTop(colPair.first, scalar1);
    builder.appendEntryTop(colPair.second, scalar2);
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
    const auto colA = findOrCreateColumn(itA.getMonomial(), mulA, builder);
    const auto colB = findOrCreateColumn(itB.getMonomial(), mulB, builder);
    const auto cmp = ring().monomialCompare
      (builder.monomialOfCol(colA), builder.monomialOfCol(colB));
    if (cmp != LT) {
      coeff = itA.getCoefficient();
      col = colA;
      ++itA;
    }
    if (cmp != GT) {
      coeff = ring().coefficientSubtract(coeff, itB.getCoefficient());
      col = colB;
      ++itB;
    }
    MATHICGB_ASSERT(coeff < std::numeric_limits<Scalar>::max());
    if (coeff != 0)
      builder.appendEntryBottom(col, static_cast<Scalar>(coeff));
  }
}
