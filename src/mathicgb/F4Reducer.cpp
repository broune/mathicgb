#include "stdinc.h"
#include "F4Reducer.hpp"

#include "F4MatrixBuilder.hpp"
#include "F4MatrixReducer.hpp"
#include <iostream>

extern int tracingLevel;

F4Reducer::F4Reducer(
  const PolyRing& ring,
  std::unique_ptr<Reducer> fallback
):
  mFallback(std::move(fallback)),
  mRing(ring) {
}

std::unique_ptr<Poly> F4Reducer::classicReduce
(const Poly& poly, const PolyBasis& basis) {
  if (tracingLevel >= 2)
    std::cerr <<
      "F4Reducer: Using fall-back reducer for single classic reduction\n";

  std::unique_ptr<Poly> p;
  p = mFallback->classicReduce(poly, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

std::unique_ptr<Poly> F4Reducer::classicTailReduce
(const Poly& poly, const PolyBasis& basis) {
  std::unique_ptr<Poly> p;
  if (tracingLevel >= 2)
    std::cerr <<
      "F4Reducer: Using fall-back reducer for single classic tail reduction\n";

  p = mFallback->classicTailReduce(poly, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

std::unique_ptr<Poly> F4Reducer::classicReduceSPoly
(const Poly& a, const Poly& b, const PolyBasis& basis) {
  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: Reducing single S-pair.\n";

  QuadMatrix qm;
  {
    F4MatrixBuilder builder(basis);
    builder.addSPolynomialToMatrix(a, b);
    builder.buildMatrixAndClear(qm);

    // there has to be something to reduce
    MATHICGB_ASSERT(qm.bottomLeft.rowCount() > 0);
  }

  SparseMatrix reduced;
  {
    F4MatrixReducer red(mThreadCount);
    red.reduce(basis.ring(), qm, reduced);
  }

  auto p = make_unique<Poly>(&basis.ring());
  if (reduced.rowCount() > 0) {
    MATHICGB_ASSERT(reduced.rowCount() == 1);
    reduced.rowToPolynomial(0, qm.rightColumnMonomials, *p);
  }
  return p;
}

void F4Reducer::classicReduceSPolySet
(std::vector<std::pair<size_t, size_t> >& spairs,
 const PolyBasis& basis,
 std::vector<std::unique_ptr<Poly> >& reducedOut)
{
  reducedOut.clear();
  if (spairs.empty())
    return;

  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: Reducing " << spairs.size() << " S-pair.\n";

  SparseMatrix reduced;
  std::vector<monomial> monomials;
  {
    QuadMatrix qm;
    {
      F4MatrixBuilder builder(basis);
      for (auto it = spairs.begin(); it != spairs.end(); ++it) {
        builder.addSPolynomialToMatrix
          (basis.poly(it->first), basis.poly(it->second));
      }
      builder.buildMatrixAndClear(qm);

      // there has to be something to reduce
      MATHICGB_ASSERT(qm.bottomLeft.rowCount() > 0);
    }
    F4MatrixReducer(mThreadCount).reduce(basis.ring(), qm, reduced);
    monomials = std::move(qm.rightColumnMonomials);
  }

  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: Extracted " << reduced.rowCount()
              << " non-zero rows\n";

  for (SparseMatrix::RowIndex row = 0; row < reduced.rowCount(); ++row) {
    auto p = make_unique<Poly>(&basis.ring());
    reduced.rowToPolynomial(row, monomials, *p);
    
    reducedOut.push_back(std::move(p));
  }
}

void F4Reducer::classicReducePolySet
(const std::vector<std::unique_ptr<Poly> >& polys,
 const PolyBasis& basis,
 std::vector<std::unique_ptr<Poly> >& reducedOut)
{
  for (auto it = polys.begin(); it != polys.end(); ++it) {
    auto reducedPoly = classicReduce(**it, basis);
    if (!reducedPoly->isZero())
      reducedOut.push_back(std::move(reducedPoly));
  }
}

Poly* F4Reducer::regularReduce
(const_monomial sig,
 const_monomial multiple,
 size_t basisElement,
 const GroebnerBasis& basis) {
  Poly* p = mFallback->regularReduce(sig, multiple, basisElement, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

void F4Reducer::setThreadCount(size_t threadCount) {
  mThreadCount = threadCount;
}

std::string F4Reducer::description() const {
  return "F4 reducer";
}

size_t F4Reducer::getMemoryUse() const {
  return 0; // @todo: implement
}
