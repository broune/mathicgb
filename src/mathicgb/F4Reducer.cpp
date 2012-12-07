#include "stdinc.h"
#include "F4Reducer.hpp"

#include "F4MatrixBuilder.hpp"
#include "F4MatrixReducer.hpp"
#include "QuadMatrix.hpp"
#include <iostream>
#include <limits>

F4Reducer::F4Reducer(const PolyRing& ring):
  mFallback(Reducer::makeReducer(Reducer::Reducer_BjarkeGeo, ring)),
  mRing(ring),
  mMemoryQuantum(0),
  mStoreToFile(""),
  mMinEntryCountForStore(0),
  mMatrixSaveCount(0) {
}

void F4Reducer::writeMatricesTo(std::string file, size_t minEntries) {
  mStoreToFile = std::move(file);
  mMinEntryCountForStore = minEntries;
  mMatrixSaveCount = 0;
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

std::unique_ptr<Poly> F4Reducer::classicReduceSPoly(
  const Poly& a,
  const Poly& b,
  const PolyBasis& basis
) {
  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: "
      "Using fall-back reducer for single classic S-pair reduction\n";
  return mFallback->classicReduceSPoly(a, b, basis);
}

void F4Reducer::classicReduceSPolySet(
  std::vector<std::pair<size_t, size_t> >& spairs,
  const PolyBasis& basis,
  std::vector<std::unique_ptr<Poly> >& reducedOut
) {
  if (spairs.size() <= 1) {
    if (tracingLevel >= 2)
      std::cerr << "F4Reducer: Using fall-back reducer for "
        << spairs.size() << " S-pairs.\n";
    mFallback->classicReduceSPolySet(spairs, basis, reducedOut);
    return;
  }
  reducedOut.clear();

  MATHICGB_ASSERT(!spairs.empty());
  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: Reducing " << spairs.size() << " S-polynomials.\n";

  SparseMatrix reduced;
  std::vector<monomial> monomials;
  {
    QuadMatrix qm;
    {
      F4MatrixBuilder builder(basis, mMemoryQuantum);
      for (auto it = spairs.begin(); it != spairs.end(); ++it) {
        builder.addSPolynomialToMatrix
          (basis.poly(it->first), basis.poly(it->second));
      }
      builder.buildMatrixAndClear(qm);

      // there has to be something to reduce
      MATHICGB_ASSERT(qm.bottomLeft.rowCount() > 0);
    }
    saveMatrix(qm);
    reduced = F4MatrixReducer(basis.ring().charac()).
      reducedRowEchelonFormBottomRight(qm);
    monomials = std::move(qm.rightColumnMonomials);
    const auto end = qm.leftColumnMonomials.end();
    for (auto it = qm.leftColumnMonomials.begin(); it != end; ++it)
      mRing.freeMonomial(*it);
  }

  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: Extracted " << reduced.rowCount()
              << " non-zero rows\n";

  for (SparseMatrix::RowIndex row = 0; row < reduced.rowCount(); ++row) {
    auto p = make_unique<Poly>(basis.ring());
    reduced.rowToPolynomial(row, monomials, *p);
    reducedOut.push_back(std::move(p));
  }
  const auto end = monomials.end();
  for (auto it = monomials.begin(); it != end; ++it)
    mRing.freeMonomial(*it);
}

void F4Reducer::classicReducePolySet
(const std::vector<std::unique_ptr<Poly> >& polys,
 const PolyBasis& basis,
 std::vector<std::unique_ptr<Poly> >& reducedOut)
{
  if (polys.size() <= 1) {
    if (tracingLevel >= 2)
      std::cerr << "F4Reducer: Using fall-back reducer for "
                << polys.size() << " polynomials.\n";
    mFallback->classicReducePolySet(polys, basis, reducedOut);
    return;
  }

  reducedOut.clear();

  MATHICGB_ASSERT(!polys.empty());
  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: Reducing " << polys.size() << " polynomials.\n";

  SparseMatrix reduced;
  std::vector<monomial> monomials;
  {
    QuadMatrix qm;
    {
      F4MatrixBuilder builder(basis, mMemoryQuantum);
      for (auto it = polys.begin(); it != polys.end(); ++it)
        builder.addPolynomialToMatrix(**it);
      builder.buildMatrixAndClear(qm);

      // there has to be something to reduce
      MATHICGB_ASSERT(qm.bottomLeft.rowCount() > 0);
    }
    saveMatrix(qm);
    reduced = F4MatrixReducer(basis.ring().charac()).
      reducedRowEchelonFormBottomRight(qm);
    monomials = std::move(qm.rightColumnMonomials);
    for (auto it = qm.leftColumnMonomials.begin();
      it != qm.leftColumnMonomials.end(); ++it)
      mRing.freeMonomial(*it);
  }

  if (tracingLevel >= 2)
    std::cerr << "F4Reducer: Extracted " << reduced.rowCount()
              << " non-zero rows\n";

  for (SparseMatrix::RowIndex row = 0; row < reduced.rowCount(); ++row) {
    auto p = make_unique<Poly>(basis.ring());
    reduced.rowToPolynomial(row, monomials, *p);
    reducedOut.push_back(std::move(p));
  }
  const auto end = monomials.end();
  for (auto it = monomials.begin(); it != end; ++it)
    mRing.freeMonomial(*it);
}

Poly* F4Reducer::regularReduce(
  const_monomial sig,
  const_monomial multiple,
  size_t basisElement,
  const GroebnerBasis& basis
) {
  if (tracingLevel >= 2)
    std::cerr <<
      "F4Reducer: Using fall-back reducer for single regular reduction\n";
  Poly* p = mFallback->regularReduce(sig, multiple, basisElement, basis);
  mSigStats = mFallback->sigStats();
  mClassicStats = mFallback->classicStats();
  return p;
}

void F4Reducer::setMemoryQuantum(size_t quantum) {
  mMemoryQuantum = quantum;
}

std::string F4Reducer::description() const {
  return "F4 reducer";
}

size_t F4Reducer::getMemoryUse() const {
  return 0; // @todo: implement
}

void F4Reducer::saveMatrix(const QuadMatrix& matrix) {
  if (mStoreToFile.empty())
    return;
  const auto entryCount = matrix.entryCount();
  if (mMinEntryCountForStore > entryCount)
    return;
  ++mMatrixSaveCount;
  std::ostringstream fileName;
  fileName << mStoreToFile << '-' << mMatrixSaveCount << ".qmat";
  if (tracingLevel > 2)
    std::cerr << "F4Reducer: Saving matrix to " << fileName.str() << '\n';
  FILE* file = fopen(fileName.str().c_str(), "wb");
  // @todo: fix leak of file on exception
  matrix.write(static_cast<SparseMatrix::Scalar>(mRing.charac()), file);
  fclose(file);
}
