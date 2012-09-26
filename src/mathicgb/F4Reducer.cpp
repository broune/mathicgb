#include "stdinc.h"
#include "F4Reducer.hpp"

#include "F4MatrixBuilder.hpp"

F4Reducer::F4Reducer(const PolyRing& ring, std::auto_ptr<Reducer> fallback):
  mFallback(fallback), mRing(ring) {
}

std::auto_ptr<Poly> F4Reducer::classicReduce
(const Poly& poly, const PolyBasis& basis) {
  std::auto_ptr<Poly> p;
  p = mFallback->classicReduce(poly, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

std::auto_ptr<Poly> F4Reducer::classicTailReduce
(const Poly& poly, const PolyBasis& basis) {
  std::auto_ptr<Poly> p;
  p = mFallback->classicTailReduce(poly, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

std::auto_ptr<Poly> F4Reducer::classicReduceSPoly
(const Poly& a, const Poly& b, const PolyBasis& basis) {
  std::auto_ptr<Poly> p;
  {
    p = mFallback->classicReduceSPoly(a, b, basis);
    mSigStats = mFallback->sigStats();
    mClassicStats = mFallback->classicStats();
  }

  SparseMatrix topLeft;
  SparseMatrix topRight;
  SparseMatrix bottomLeft;
  SparseMatrix bottomRight;
  std::vector<monomial> monomialsOfRightColumns;
  {
    F4MatrixBuilder builder(basis);
    builder.addTwoRowsForSPairToMatrix(a, b);
    builder.buildMatricesAndClear
      (topLeft, topRight, bottomLeft, bottomRight, monomialsOfRightColumns);
  }

  return p;
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
  return 0;
}
