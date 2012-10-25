#include "stdinc.h"
#include "F4MatrixBuilder.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

F4MatrixBuilder::F4MatrixBuilder(const PolyBasis& basis, const int threadCount):
  mThreadCount(threadCount),
  mBasis(basis),
  mBuilder(basis.ring()),
  mTmp(basis.ring().allocMonomial())
{
  MATHICGB_ASSERT(threadCount >= 1);
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
        auto col = mBuilder.findColumn(mono);
        if (!col.valid())
          col = this->createColumn(mono, mBuilder);
        mBuilder.appendEntryBottom(col,
          static_cast<QuadMatrixBuilder::Scalar>(c));
	  }
    }
    ring().freeMonomial(monoA);
    ring().freeMonomial(monoB);
    mBuilder.rowDoneBottomLeftAndRight();
  }

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
      const monomial tmp = td.tmp;
#else
      const monomial tmp = mTmp;
      QuadMatrixBuilder& builder = mBuilder;
#endif
      const RowTask task = currentTasks[taskIndex];
      MATHICGB_ASSERT(ring().hashValid(task.multiple));

      appendRowTop(task.multiple, *task.poly, builder, tmp);
    }
  }

#ifdef _OPENMP
#pragma omp flush
  MATHICGB_ASSERT(mThreadCount > 1 || threadData.back().builder == &mBuilder);
  if (mThreadCount > 1) {
    for (auto it = threadData.begin(); it != threadData.end(); ++it) {
      MATHICGB_ASSERT(&mBuilder != it->builder);
      QuadMatrix qm;
      it->builder->buildMatrixAndClear(qm);
      mBuilder.takeRowsFrom(std::move(qm));
      delete it->builder;
      delete[] it->tmp.unsafeGetRepresentation();
    }
  }
  threadData.clear();
#endif


  mBuilder.sortColumnsLeftRightParallel(mBasis.order(), mThreadCount);
  mBuilder.buildMatrixAndClear(matrix);
}

F4MatrixBuilder::LeftRightColIndex
F4MatrixBuilder::createColumn(
  const const_monomial mono,
  QuadMatrixBuilder& builder
) {
    QuadMatrixBuilder::LeftRightColIndex newCol;

#pragma omp critical
  {
    MATHICGB_ASSERT(ring().hashValid(mono));
    MATHICGB_ASSERT(!builder.findColumn(mono).valid());
#ifndef _OPENMP
    MATHICGB_ASSERT(&builder == &mBuilder);
#endif

#ifdef _OPENMP
    // builder does not have a column for mono but mBuilder might
    if (&mBuilder != &builder)
      newCol = mBuilder.findColumn(mono);
    if (!newCol.valid())
#endif
    {
      // we need to add a new column to mBuilder

      // look for a reducer of mono
      size_t reducerIndex = mBasis.divisor(mono);
      if (reducerIndex == static_cast<size_t>(-1)) {
        newCol = mBuilder.createColumnRight(mono);
      } else {
        newCol = mBuilder.createColumnLeft(mono);

        // schedule the reducer to be added as a row
        RowTask task;
        task.poly = &mBasis.poly(reducerIndex);
        task.multiple = ring().allocMonomial();
        ring().monomialDivideToNegative
          (mono, task.poly->getLeadMonomial(), task.multiple);
        MATHICGB_ASSERT(ring().hashValid(task.multiple));
        mTodo.push_back(task);
      }
    }
    MATHICGB_ASSERT(newCol.valid());
    MATHICGB_ASSERT(ring().monomialEQ(mBuilder.monomialOfCol(newCol), mono));
    MATHICGB_ASSERT(newCol == mBuilder.findColumn(mono));

#ifdef _OPENMP
    // If we are running in parallel we need to catch builder up to date with
    // with respect to the columns in mBuilder. We have just ensured that
    // mBuilder has a column for mono so now builder will also get that column.

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

      MATHICGB_ASSERT(ring().monomialEQ(mBuilder.monomialOfCol(mBuilder.findColumn(mono)), mono));
      MATHICGB_ASSERT(newCol == mBuilder.findColumn(mono));

      MATHICGB_ASSERT(builder.findColumn(mono).valid());
      MATHICGB_ASSERT(ring().monomialEQ(builder.monomialOfCol(builder.findColumn(mono)), mono));
      MATHICGB_ASSERT(newCol == builder.findColumn(mono));
    }
#endif
  }
  return newCol;
}

void F4MatrixBuilder::appendRowTop(const const_monomial multiple, const Poly& poly, QuadMatrixBuilder& builder, monomial tmp) {
  Poly::const_iterator end = poly.end();
  for (Poly::const_iterator it = poly.begin(); it != end; ++it) {
	MATHICGB_ASSERT(it.getCoefficient() <
      std::numeric_limits<QuadMatrixBuilder::Scalar>::max());
    auto col = builder.findColumnProduct(it.getMonomial(), multiple);
    if (!col.valid()) {
      ring().monomialMult(it.getMonomial(), multiple, tmp);
      col = createColumn(tmp, builder);
    }
    const auto scalar =
      static_cast<QuadMatrixBuilder::Scalar>(it.getCoefficient());
    builder.appendEntryTop(col, scalar);
  }
  builder.rowDoneTopLeftAndRight();
}
