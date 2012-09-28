#include "stdinc.h"
#include "F4Reducer.hpp"

#include "F4MatrixBuilder.hpp"
#include "F4MatrixReducer.hpp"

F4Reducer::F4Reducer(const PolyRing& ring, std::unique_ptr<Reducer> fallback):
  mFallback(std::move(fallback)), mRing(ring) {
}

std::unique_ptr<Poly> F4Reducer::classicReduce
(const Poly& poly, const PolyBasis& basis) {
  std::unique_ptr<Poly> p;
  p = mFallback->classicReduce(poly, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

std::unique_ptr<Poly> F4Reducer::classicTailReduce
(const Poly& poly, const PolyBasis& basis) {
  std::unique_ptr<Poly> p;
  p = mFallback->classicTailReduce(poly, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

std::unique_ptr<Poly> F4Reducer::classicReduceSPoly
(const Poly& a, const Poly& b, const PolyBasis& basis) {
  QuadMatrix qm;
  {
    F4MatrixBuilder builder(basis);
    builder.addTwoRowsForSPairToMatrix(a, b);
    builder.buildMatrixAndClear(qm);

    // there has to be something to reduce
    MATHICGB_ASSERT(qm.bottomLeft.rowCount() > 0);
  }

  SparseMatrix reduced;
  {
    F4MatrixReducer red;
    red.reduce(basis.ring(), qm, reduced);
  }

  auto p = make_unique<Poly>(&basis.ring());
  if (reduced.rowCount() > 0) {
    MATHICGB_ASSERT(reduced.rowCount() == 1);
    reduced.rowToPolynomial(0, qm.rightColumnMonomials, *p);
  }
  return p;
}

void F4Reducer::classicReduceSPolyGroup
(std::vector<std::pair<size_t, size_t> >& spairs,
 const PolyBasis& basis,
 std::vector<std::unique_ptr<Poly> >& reducedOut)
{
  reducedOut.clear();
  if (spairs.empty())
    return;

  SparseMatrix reduced;
  std::vector<monomial> monomials;
  {
    QuadMatrix qm;
    {
      F4MatrixBuilder builder(basis);
      for (auto it = spairs.begin(); it != spairs.end(); ++it) {
        builder.addTwoRowsForSPairToMatrix
          (basis.poly(it->first), basis.poly(it->second));
      }
      builder.buildMatrixAndClear(qm);

      // there has to be something to reduce
      MATHICGB_ASSERT(qm.bottomLeft.rowCount() > 0);
    }
    F4MatrixReducer().reduce(basis.ring(), qm, reduced);
    monomials = std::move(qm.rightColumnMonomials);
  }

  for (SparseMatrix::RowIndex row = 0; row < reduced.rowCount(); ++row) {
    auto p = make_unique<Poly>(&basis.ring());
    reduced.rowToPolynomial(row, monomials, *p);
    reducedOut.push_back(std::move(p));
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

std::string F4Reducer::description() const {
  return "F4 reducer";
}

size_t F4Reducer::getMemoryUse() const {
  return 0; // @todo: implement
}
