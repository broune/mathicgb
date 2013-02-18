#include "stdinc.h"
#include "F4MatrixBuilder2.hpp"

#include "LogDomain.hpp"
#include "F4MatrixProjection.hpp"
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
  mMap(basis.ring())
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
  task.poly = &polyA;
  task.sPairPoly = &polyB;
  mTodo.push_back(task);
}

void F4MatrixBuilder2::addPolynomialToMatrix(const Poly& poly) {
  if (poly.isZero())
    return;

  RowTask task = {};
  task.poly = &poly;
  mTodo.push_back(task);
}

void F4MatrixBuilder2::addPolynomialToMatrix
(const_monomial multiple, const Poly& poly) {
  MATHICGB_ASSERT(ring().hashValid(multiple));
  if (poly.isZero())
    return;

  RowTask task = {};
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

  // Process pending rows until we are done. Note that the methods
  // we are calling here can add more pending items.

  struct ThreadData {
    F4ProtoMatrix block;
    monomial tmp1;
    monomial tmp2;
  };

  tbb::enumerable_thread_specific<ThreadData> threadData([&](){  
    ThreadData data;
    {
      tbb::mutex::scoped_lock guard(mCreateColumnLock);
      data.tmp1 = ring().allocMonomial();
      data.tmp2 = ring().allocMonomial();
    }
    return std::move(data);
  });

  // Construct the matrix as pre-blocks
  tbb::parallel_do(mTodo.begin(), mTodo.end(),
    [&](const RowTask& task, TaskFeeder& feeder)
  {
    auto& data = threadData.local();
    const auto& poly = *task.poly;

    if (task.sPairPoly != 0) {
      ring().monomialColons(
        poly.getLeadMonomial(),
        task.sPairPoly->getLeadMonomial(),
        data.tmp2,
        data.tmp1
      );
      appendRowSPair
        (&poly, data.tmp1, task.sPairPoly, data.tmp2, data.block, feeder);
      return;
    }
    if (task.desiredLead.isNull())
      ring().monomialSetIdentity(data.tmp1);
    else
      ring().monomialDivide
        (task.desiredLead, poly.getLeadMonomial(), data.tmp1);
    MATHICGB_ASSERT(ring().hashValid(data.tmp1));
    appendRow(data.tmp1, *task.poly, data.block, feeder);
  });
  MATHICGB_ASSERT(!threadData.empty()); // as mTodo empty causes early return

  // Free the monomials from all the tasks
  const auto todoEnd = mTodo.end();
  for (auto it = mTodo.begin(); it != todoEnd; ++it)
    if (!it->desiredLead.isNull())
      ring().freeMonomial(it->desiredLead);
  mTodo.clear();

  // Move the proto-matrices across all threads into the projection.
  F4MatrixProjection projection
    (ring(), static_cast<ColIndex>(mMap.entryCount()));
  const auto end = threadData.end();
  for (auto it = threadData.begin(); it != end; ++it) {
    projection.addProtoMatrix(std::move(it->block));
    ring().freeMonomial(it->tmp1);
    ring().freeMonomial(it->tmp2);
  }

  // Sort columns by monomial and tell the projection of the resulting order
  MonomialMap<ColIndex>::Reader reader(mMap);
  typedef std::pair<ColIndex, ConstMonomial> IndexMono;
  std::vector<IndexMono> columns(reader.begin(), reader.end());
  const auto cmp = [&](const IndexMono& a, const IndexMono b) {
    return ring().monomialLT(b.second, a.second);
  };
  tbb::parallel_sort(columns.begin(), columns.end(), cmp);

  const auto colEnd = columns.end();
  for (auto it = columns.begin(); it != colEnd; ++it) {
    const auto p = *it;
    projection.addColumn(p.first, p.second, mIsColumnToLeft[p.first]);
  }

  quadMatrix = projection.makeProjectionAndClear();
  threadData.clear();

#ifdef MATHICGB_DEBUG
  for (size_t side = 0; side < 2; ++side) {
    auto& monos = side == 0 ?
      quadMatrix.leftColumnMonomials : quadMatrix.rightColumnMonomials;
    for (auto it = monos.begin(); it != monos.end(); ++it) {
      MATHICGB_ASSERT(!it->isNull());
    }
  }
  for (RowIndex row = 0; row < quadMatrix.topLeft.rowCount(); ++row) {
    MATHICGB_ASSERT(quadMatrix.topLeft.entryCountInRow(row) > 0);
    MATHICGB_ASSERT(quadMatrix.topLeft.leadCol(row) == row);
  }
  MATHICGB_ASSERT(quadMatrix.debugAssertValid());
#endif
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
  const size_t reducerIndex = mBasis.classicReducer(mTmp);
  const bool insertLeft = (reducerIndex != static_cast<size_t>(-1));

  // Create the new left or right column
  if (mIsColumnToLeft.size() >= std::numeric_limits<ColIndex>::max())
    throw std::overflow_error("Too many columns in QuadMatrix");
  const auto newIndex = static_cast<ColIndex>(mIsColumnToLeft.size());
  const auto inserted = mMap.insert(std::make_pair(mTmp, newIndex));
  mIsColumnToLeft.push_back(insertLeft);

  // schedule new task if we found a reducer
  if (insertLeft) {
    RowTask task = {};
    task.poly = &mBasis.poly(reducerIndex);
    task.desiredLead = inserted.first.second.castAwayConst();
    feeder.add(task);
  }

  return std::make_pair(*inserted.first.first, inserted.first.second);
}

void F4MatrixBuilder2::appendRow(
  const const_monomial multiple,
  const Poly& poly,
  F4ProtoMatrix& block,
  TaskFeeder& feeder
) {
  MATHICGB_ASSERT(!multiple.isNull());

  const auto begin = poly.begin();
  const auto end = poly.end();
  const auto count = std::distance(begin, end);
  MATHICGB_ASSERT(count < std::numeric_limits<ColIndex>::max());
  auto indices = block.makeRowWithTheseScalars(poly);

  auto it = begin;
  if ((count % 2) == 1) {
    ColReader reader(mMap);
    const auto col = findOrCreateColumn
      (it.getMonomial(), multiple, reader, feeder);
	MATHICGB_ASSERT(it.getCoefficient() < std::numeric_limits<Scalar>::max());
    MATHICGB_ASSERT(it.getCoefficient());
    //matrix.appendEntry(col.first, static_cast<Scalar>(it.getCoefficient()));
    *indices = col.first;
    ++indices;
    //*row.first++ = col.first;
    //*row.second++ = static_cast<Scalar>(it.getCoefficient());
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

    //matrix.appendEntry(*colPair.first, scalar1);
    //matrix.appendEntry(*colPair.second, scalar2);

    *indices = *colPair.first;
    ++indices;
    *indices = *colPair.second;
    ++indices;

    //*row.first++ = *colPair.first;
    //*row.second++ = scalar1;
    //*row.first++ = *colPair.second;
    //*row.second++ = scalar2;

    it = ++it2;
  }
  //matrix.rowDone();
}

void F4MatrixBuilder2::appendRowSPair(
  const Poly* poly,
  monomial multiply,
  const Poly* sPairPoly,
  monomial sPairMultiply,
  F4ProtoMatrix& block,
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

  // @todo: handle overflow of termCount addition here
  MATHICGB_ASSERT(poly->termCount() + sPairPoly->termCount() - 2 <=
    std::numeric_limits<ColIndex>::max());
  const auto maxCols =
    static_cast<ColIndex>(poly->termCount() + sPairPoly->termCount() - 2);
  auto row = block.makeRow(maxCols);
  const auto indicesBegin = row.first;

  const ColReader colMap(mMap);

  const const_monomial mulA = multiply;
  const const_monomial mulB = sPairMultiply;
  while (itB != endB && itA != endA) {
    const auto colA = findOrCreateColumn
      (itA.getMonomial(), mulA, colMap, feeder);
    const auto colB = findOrCreateColumn
      (itB.getMonomial(), mulB, colMap, feeder);
    const auto cmp = ring().monomialCompare(colA.second, colB.second);

    coefficient coeff = 0;
    ColIndex col;
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
      //matrix.appendEntry(col, static_cast<Scalar>(coeff));
      *row.first++ = col;
      *row.second++ = static_cast<Scalar>(coeff);
    }
  }

  for (; itA != endA; ++itA) {
    const auto colA = findOrCreateColumn
      (itA.getMonomial(), mulA, colMap, feeder);
    //matrix.appendEntry(colA.first, static_cast<Scalar>(itA.getCoefficient()));
    *row.first++ = colA.first;
    *row.second++ = static_cast<Scalar>(itA.getCoefficient());
  }

  for (; itB != endB; ++itB) {
    const auto colB = findOrCreateColumn
      (itB.getMonomial(), mulB, colMap, feeder);
    const auto negative = ring().coefficientNegate(itB.getCoefficient());
    //matrix.appendEntry(colB.first, static_cast<Scalar>(negative));
    *row.first++ = colB.first;
    *row.second++ = static_cast<Scalar>(negative);
  }

  const auto toRemove =
    maxCols - static_cast<ColIndex>(row.first - indicesBegin);
  block.removeLastEntries(block.rowCount() - 1, toRemove);
}
