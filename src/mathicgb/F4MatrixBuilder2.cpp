#include "stdinc.h"
#include "F4MatrixBuilder2.hpp"

#include "LogDomain.hpp"
#include <tbb/tbb.h>

MATHICGB_DEFINE_LOG_DOMAIN(
  F4MatrixBuild2,
  "Displays statistics about F4 matrix construction."
);

MATHICGB_NO_INLINE
std::pair<F4MatrixBuilder2::ColIndex, ConstMonomial>
F4MatrixBuilder2::findOrCreateColumn(
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
std::pair<F4MatrixBuilder2::ColIndex, ConstMonomial>
F4MatrixBuilder2::findOrCreateColumn(
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
void F4MatrixBuilder2::createTwoColumns(
  const const_monomial monoA1,
  const const_monomial monoA2,
  const const_monomial monoB,
  TaskFeeder& feeder
) {
  createColumn(monoA1, monoB, feeder);
  createColumn(monoA2, monoB, feeder);
}

F4MatrixBuilder2::F4MatrixBuilder2(
  const PolyBasis& basis,
  const size_t memoryQuantum
):
  mMemoryQuantum(memoryQuantum),
  mTmp(basis.ring().allocMonomial()),
  mBasis(basis),
  mMap(basis.ring()),
  mLeftColCount(0),
  mRightColCount(0)
{
  // This assert has to be _NO_ASSUME since otherwise the compiler will assume
  // that the error checking branch here cannot be taken and optimize it away.
  const Scalar maxScalar = std::numeric_limits<Scalar>::max();
  MATHICGB_ASSERT_NO_ASSUME(ring().charac() <= maxScalar);
  if (ring().charac() > maxScalar)
    mathic::reportInternalError("F4MatrixBuilder2: too large characteristic.");
}

void F4MatrixBuilder2::addSPolynomialToMatrix(
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

void F4MatrixBuilder2::addPolynomialToMatrix(const Poly& poly) {
  if (poly.isZero())
    return;

  RowTask task = {};
  task.addToTop = false;
  task.poly = &poly;
  mTodo.push_back(task);
}

void F4MatrixBuilder2::addPolynomialToMatrix
(const_monomial multiple, const Poly& poly) {
  MATHICGB_ASSERT(ring().hashValid(multiple));
  if (poly.isZero())
    return;

  RowTask task = {};
  task.addToTop = false;
  task.poly = &poly;
  task.desiredLead = ring().allocMonomial();
  ring().monomialMult(poly.getLeadMonomial(), multiple, task.desiredLead);
  MATHICGB_ASSERT(ring().hashValid(task.desiredLead));

  MATHICGB_ASSERT(task.sPairPoly == 0);
  mTodo.push_back(task);
}

void F4MatrixBuilder2::buildMatrixAndClear(QuadMatrix& quadMatrix) {
  MATHICGB_LOG_TIME(F4MatrixBuild2) <<
    "\n***** Constructing matrix *****\n";

  if (mTodo.empty()) {
    quadMatrix = QuadMatrix();
    quadMatrix.ring = &ring();
    return;
  }

  // todo: prefer sparse/old reducers among the inputs.

  // Process pending rows until we are done. Note that the methods
  // we are calling here can add more pending items.

  struct ThreadData {
    SparseMatrix topMatrix;
    SparseMatrix bottomMatrix;
    monomial tmp1;
    monomial tmp2;
  };

  tbb::enumerable_thread_specific<ThreadData> threadData([&](){  
    ThreadData data = {
      SparseMatrix(mMemoryQuantum),
      SparseMatrix(mMemoryQuantum)
    };
    {
      tbb::mutex::scoped_lock guard(mCreateColumnLock);
      data.tmp1 = ring().allocMonomial();
      data.tmp2 = ring().allocMonomial();
    }
    return std::move(data);
  });

  tbb::parallel_do(mTodo.begin(), mTodo.end(),
    [&](const RowTask& task, TaskFeeder& feeder)
  {
    auto& data = threadData.local();
    const auto& poly = *task.poly;

    if (task.sPairPoly != 0) {
      MATHICGB_ASSERT(!task.addToTop);
      ring().monomialColons(
        poly.getLeadMonomial(),
        task.sPairPoly->getLeadMonomial(),
        data.tmp2,
        data.tmp1
      );
      appendRowBottom
        (&poly, data.tmp1, task.sPairPoly, data.tmp2, data.bottomMatrix, feeder);
      MATHICGB_ASSERT(data.bottomMatrix.rowCount() == 0 || !data.bottomMatrix.emptyRow(data.bottomMatrix.rowCount() - 1));
      return;
    }
    if (task.desiredLead.isNull())
      ring().monomialSetIdentity(data.tmp1);
    else
      ring().monomialDivide
        (task.desiredLead, poly.getLeadMonomial(), data.tmp1);
    MATHICGB_ASSERT(ring().hashValid(data.tmp1));
    if (task.addToTop) {
      appendRowTop(data.tmp1, *task.poly, data.topMatrix, feeder);
      MATHICGB_ASSERT(!data.topMatrix.emptyRow(data.topMatrix.rowCount() - 1));
    } else {
      appendRowBottom
        (data.tmp1, false, poly.begin(), poly.end(), data.bottomMatrix, feeder);
      MATHICGB_ASSERT(data.bottomMatrix.rowCount() == 0 || !data.bottomMatrix.emptyRow(data.bottomMatrix.rowCount() - 1));
    }
  });
  MATHICGB_ASSERT(!threadData.empty()); // as mTodo empty causes early return

  auto topMatrix = std::move(threadData.begin()->topMatrix);
  auto bottomMatrix = std::move(threadData.begin()->bottomMatrix);
  const auto end = threadData.end();
  for (auto it = threadData.begin(); it != end; ++it) {
    if (it != threadData.begin()) {
      topMatrix.takeRowsFrom(std::move(it->topMatrix));
      bottomMatrix.takeRowsFrom(std::move(it->bottomMatrix));
    }
    ring().freeMonomial(it->tmp1);
    ring().freeMonomial(it->tmp2);
  }
  threadData.clear();

  MATHICGB_ASSERT(bottomMatrix.rowCount() == mTodo.size());

  // Free the monomials from all the tasks
  const auto todoEnd = mTodo.end();
  for (auto it = mTodo.begin(); it != todoEnd; ++it)
    if (!it->desiredLead.isNull())
      ring().freeMonomial(it->desiredLead);
  mTodo.clear();

  {
    ColReader reader(mMap);
    quadMatrix.leftColumnMonomials.swap(Monomials(mLeftColCount));
    quadMatrix.rightColumnMonomials.swap(Monomials(mRightColCount));
    const auto end = reader.end();
    for (auto it = reader.begin(); it != end; ++it) {
      const auto p = *it;
      monomial copy = ring().allocMonomial();
      ring().monomialCopy(p.second, copy);
      const auto index = p.first;
      auto translated = mTranslate[index];
      auto& monos = translated.left ?
        quadMatrix.leftColumnMonomials : quadMatrix.rightColumnMonomials;
      MATHICGB_ASSERT(translated.index < monos.size());
      MATHICGB_ASSERT(monos[translated.index].isNull());
      monos[translated.index] = copy;
    }
  }

  topMatrix.takeRowsFrom(std::move(bottomMatrix));
  auto& matrix = topMatrix;

#ifdef MATHICGB_DEBUG
  for (size_t row = 0; row < matrix.rowCount(); ++row) {
    MATHICGB_ASSERT(!matrix.emptyRow(row));
  }
#endif

  // Decide which rows are reducers (top) and which are reducees (bottom).
  const auto rowCount = matrix.rowCount();
  std::vector<RowIndex> reducerRows(mLeftColCount, rowCount);
  std::vector<RowIndex> reduceeRows;
  for (RowIndex row = 0; row < rowCount; ++row) {
    MATHICGB_ASSERT(!matrix.emptyRow(row));
    // Determine leading (minimum index) left entry.
    const auto lead = [&] {
      MATHICGB_ASSERT(mTranslate.size() <= std::numeric_limits<ColIndex>::max());
      const auto end = matrix.rowEnd(row);
      for (auto it = matrix.rowBegin(row); it != end; ++it) {
        MATHICGB_ASSERT(it.index() < mTranslate.size());
        if (mTranslate[it.index()].left)
          return mTranslate[it.index()].index;
      }
      return mLeftColCount; // no left entries at all
    }();
    if (!mTranslate[matrix.leadCol(row)].left) {
      reduceeRows.push_back(row); // no left entries
      continue;
    }

    // decide if this is a reducer or reducee row
    if (lead == mLeftColCount) {
      reduceeRows.push_back(row); // no left entries
      continue;
    }
    auto& reducer = reducerRows[lead];
    if (reducer == rowCount)
      reducer = row; // row is first reducer with this lead
    else if (matrix.entryCountInRow(reducer) > matrix.entryCountInRow(row)) {
      reduceeRows.push_back(reducer); // row sparser (=better) reducer
      reducer = row;
    } else
      reduceeRows.push_back(row);
  }

  MATHICGB_ASSERT(reducerRows.size() == mLeftColCount);
  MATHICGB_ASSERT(reduceeRows.size() == matrix.rowCount() - mLeftColCount);
#ifdef MATHICGB_DEBUG
  for (size_t  i = 0; i < reducerRows.size(); ++i) {
    const auto row = reducerRows[i];
    MATHICGB_ASSERT(row < matrix.rowCount());
    MATHICGB_ASSERT(!matrix.emptyRow(row));
    MATHICGB_ASSERT(mTranslate[matrix.leadCol(row)].left);
    MATHICGB_ASSERT(mTranslate[matrix.leadCol(row)].index == i);
  }
  std::vector<RowIndex> tmp;
  tmp.insert(tmp.end(), reducerRows.begin(), reducerRows.end());
  tmp.insert(tmp.end(), reduceeRows.begin(), reduceeRows.end());
  std::sort(tmp.begin(), tmp.end());
  MATHICGB_ASSERT(tmp.size() == matrix.rowCount());
  for (size_t i = 0; i < tmp.size(); ++i) {
    MATHICGB_ASSERT(tmp[i] == i);
  }
#endif
  
  quadMatrix.ring = &ring();
  auto splitLeftRight = [this](
    SparseMatrix& from,
    const std::vector<RowIndex>& fromRows,
    const bool makeLeftUnitary,
    SparseMatrix& left,
    SparseMatrix& right
  ) {
    left.clear();
    right.clear();
    const auto fromRowsEnd = fromRows.end();
    for (
      auto fromRowsIt = fromRows.begin();
      fromRowsIt != fromRowsEnd;
      ++fromRowsIt
    ) {
      const auto row = *fromRowsIt;
      const auto fromEnd = from.rowEnd(row);
      auto fromIt = from.rowBegin(row);
      MATHICGB_ASSERT(!from.emptyRow(row));
      if (
        makeLeftUnitary &&
        (!mTranslate[fromIt.index()].left || fromIt.scalar() != 1)
      ) {
        auto firstLeft = fromIt;
        while (!mTranslate[firstLeft.index()].left) {
          // We cannot make the left part unitary if the left part is a zero
          // row, so makeLeftUnitary implies no zero left rows.
          MATHICGB_ASSERT(firstLeft != fromEnd);
          ++firstLeft;
        }
        if (firstLeft.scalar() != 1) {
          const auto modulus = static_cast<Scalar>(ring().charac());
          const auto inverse = modularInverse(firstLeft.scalar(), modulus);
          for (auto it = fromIt; it != fromEnd; ++it)
            it.setScalar(modularProduct(it.scalar(), inverse, modulus));
          MATHICGB_ASSERT(firstLeft.scalar() == 1);
        }
      }
      for (; fromIt != fromEnd; ++fromIt) {
        MATHICGB_ASSERT(fromIt.index() < mTranslate.size());
        const auto translated = mTranslate[fromIt.index()];
        if (translated.left)
          left.appendEntry(translated.index, fromIt.scalar());
        else
          right.appendEntry(translated.index, fromIt.scalar());
      }
      left.rowDone();
      right.rowDone();
      MATHICGB_ASSERT(left.rowCount() == right.rowCount());
      MATHICGB_ASSERT(!makeLeftUnitary || !left.emptyRow(left.rowCount() - 1));
      MATHICGB_ASSERT
        (!makeLeftUnitary || left.rowBegin(left.rowCount() - 1).scalar() == 1);
    }
  };
  splitLeftRight
    (matrix, reducerRows, true, quadMatrix.topLeft, quadMatrix.topRight);
  splitLeftRight(
    matrix,
    reduceeRows,
    false,
    quadMatrix.bottomLeft,
    quadMatrix.bottomRight
  );

#ifdef MATHICGB_DEBUG
  for (size_t side = 0; side < 2; ++side) {
    auto& monos = side == 0 ?
      quadMatrix.leftColumnMonomials : quadMatrix.rightColumnMonomials;
    for (auto it = monos.begin(); it != monos.end(); ++it) {
      MATHICGB_ASSERT(!it->isNull());
    }
  }
  for (size_t row = 0; row < quadMatrix.topLeft.rowCount(); ++row) {
    MATHICGB_ASSERT(quadMatrix.topLeft.leadCol(row) == row);
  }
#endif

  quadMatrix.sortColumnsLeftRightParallel();
#ifdef MATHICGB_DEBUG
  MATHICGB_ASSERT(quadMatrix.debugAssertValid());
#endif

  mMap.clearNonConcurrent();
}

std::pair<F4MatrixBuilder2::ColIndex, ConstMonomial>
F4MatrixBuilder2::createColumn(
  const const_monomial monoA,
  const const_monomial monoB,
  TaskFeeder& feeder
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
  if (!ring().monomialHasAmpleCapacity(mTmp))
    mathic::reportError("Monomial exponent overflow in F4MatrixBuilder2.");
  MATHICGB_ASSERT(ring().hashValid(mTmp));

  // look for a reducer of mTmp
  const size_t reducerIndex = mBasis.divisor(mTmp);
  const bool insertLeft = (reducerIndex != static_cast<size_t>(-1));

  MATHICGB_ASSERT(mLeftColCount + mRightColCount == mTranslate.size());

  // Create the new left or right column
  auto& colCount = insertLeft ? mLeftColCount : mRightColCount;
  if (colCount == std::numeric_limits<ColIndex>::max())
    throw std::overflow_error("Too many columns in QuadMatrix");
  const auto newIndex = static_cast<ColIndex>(mTranslate.size());
  const auto inserted = mMap.insert(std::make_pair(mTmp, newIndex));
  Translated translated = {colCount, insertLeft};
  mTranslate.push_back(translated);
  ++colCount;

  MATHICGB_ASSERT(mLeftColCount + mRightColCount == mTranslate.size());

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

void F4MatrixBuilder2::appendRowBottom(
  const_monomial multiple,
  const bool negate,
  const Poly::const_iterator begin,
  const Poly::const_iterator end,
  SparseMatrix& matrix,
  TaskFeeder& feeder
) {
  // todo: eliminate the code-duplication between here and appendRowTop.
  MATHICGB_ASSERT(!multiple.isNull());
  MATHICGB_ASSERT(&matrix != 0);

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
    matrix.appendEntry(*col.first, static_cast<Scalar>(maybeNegated));
  }
  matrix.rowDone();
}

void F4MatrixBuilder2::appendRowTop(
  const const_monomial multiple,
  const Poly& poly,
  SparseMatrix& matrix,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!multiple.isNull());
  MATHICGB_ASSERT(&poly != 0);
  MATHICGB_ASSERT(&matrix != 0);

  auto it = poly.begin();
  const auto end = poly.end();
  if ((std::distance(it, end) % 2) == 1) {
    ColReader reader(mMap);
    const auto col = findOrCreateColumn
      (it.getMonomial(), multiple, reader, feeder);
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient());
    matrix.appendEntry(col.first, static_cast<Scalar>(it.getCoefficient()));
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

    matrix.appendEntry(*colPair.first, scalar1);
    matrix.appendEntry(*colPair.second, scalar2);
    it = ++it2;
  }
  matrix.rowDone();
}

void F4MatrixBuilder2::appendRowBottom(
  const Poly* poly,
  monomial multiply,
  const Poly* sPairPoly,
  monomial sPairMultiply,
  SparseMatrix& matrix,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!poly->isZero());
  MATHICGB_ASSERT(!multiply.isNull());
  MATHICGB_ASSERT(ring().hashValid(multiply));
  MATHICGB_ASSERT(sPairPoly != 0);
  Poly::const_iterator itA = poly->begin();
  const Poly::const_iterator endA = poly->end();

  MATHICGB_ASSERT(!sPairPoly->isZero());
  MATHICGB_ASSERT(!sPairMultiply.isNull());
  MATHICGB_ASSERT(ring().hashValid(sPairMultiply));
  Poly::const_iterator itB = sPairPoly->begin();
  Poly::const_iterator endB = sPairPoly->end();

  // skip leading terms since they cancel
  MATHICGB_ASSERT(itA.getCoefficient() == itB.getCoefficient());
  ++itA;
  ++itB;

  const ColReader colMap(mMap);

  const const_monomial mulA = multiply;
  const const_monomial mulB = sPairMultiply;
  bool zero = true;
  while (true) {
    // Watch out: we are depending on appendRowBottom to finish the row, so
    // if you decide not to call that function in case
    // (itA == itA && itB == endB) then you need to do that yourself.
    if (itB == endB) {
      if (!zero || itA != endA)
        appendRowBottom(mulA, false, itA, endA, matrix, feeder);
      break;
    }
    if (itA == endA) {
      appendRowBottom(mulB, true, itB, endB, matrix, feeder);
      break;
    }

    coefficient coeff = 0;
    ColIndex col;
    const auto colA = findOrCreateColumn
      (itA.getMonomial(), mulA, colMap, feeder);
    const auto colB = findOrCreateColumn
      (itB.getMonomial(), mulB, colMap, feeder);
    const auto cmp = ring().monomialCompare(colA.second, colB.second);
    //const auto cmp = ring().monomialCompare
    //  (matrix.monomialOfCol(colA), matrix.monomialOfCol(colB));
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
    if (coeff != 0) {
      zero = false;
      matrix.appendEntry(col, static_cast<Scalar>(coeff));
    }
  }
  MATHICGB_ASSERT
    (matrix.rowCount() == 0 || !matrix.emptyRow(matrix.rowCount() - 1));
}
