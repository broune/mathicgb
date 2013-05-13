// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "QuadMatrixBuilder.hpp"

#include "QuadMatrix.hpp"
#include <mathic.h>
#include <sstream>

MATHICGB_NAMESPACE_BEGIN

QuadMatrixBuilder::QuadMatrixBuilder(
  const PolyRing& ring,
  Map& map,
  MonomialsType& monomialsLeft,
  MonomialsType& monomialsRight,
  const size_t memoryQuantum
):
  mMonomialToCol(map),
  mTopLeft(memoryQuantum),
  mTopRight(memoryQuantum),
  mBottomLeft(memoryQuantum),
  mBottomRight(memoryQuantum),
  mMonomialsLeft(monomialsLeft),
  mMonomialsRight(monomialsRight)
{}

void QuadMatrixBuilder::takeRowsFrom(QuadMatrix&& matrix) {
  MATHICGB_ASSERT(&ring() == matrix.ring);
  MATHICGB_ASSERT(matrix.debugAssertValid());

  mTopLeft.takeRowsFrom(::std::move(matrix.topLeft));
  mTopRight.takeRowsFrom(::std::move(matrix.topRight));
  mBottomLeft.takeRowsFrom(::std::move(matrix.bottomLeft));
  mBottomRight.takeRowsFrom(::std::move(matrix.bottomRight));
}


namespace {
  /// Creates a column and updates the associated data structures that
  /// are passed in. Copies mono - ownership is not taken over. The
  /// purpose of this function is to avoid code duplication. It is a            
  /// template in order to avoid referring to private types of
  /// QuadMatrixBuilder.
  template<class ToMono, class ToCol>
  ::std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial>
  createCol(
    const_monomial mono,
    SparseMatrix& top,
    SparseMatrix& bottom,
    ToMono& toMonomial,
    ToCol& toCol,
    const PolyRing& ring,
    const bool left
  ) {
    MATHICGB_ASSERT(typename ToCol::Reader(toCol).find(mono).first == 0);

    const auto colCount =
      static_cast<QuadMatrixBuilder::ColIndex>(toMonomial.size());
    if (colCount == ::std::numeric_limits<QuadMatrixBuilder::ColIndex>::max())
      throw ::std::overflow_error("Too many columns in QuadMatrixBuilder");

    toMonomial.push_back(0); // allocate memory now to avoid bad_alloc later
    monomial copied = ring.allocMonomial();
    ring.monomialCopy(mono, copied);
    ::std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial> p;
    try {
      auto inserted = toCol.insert(::std::make_pair(
          copied, QuadMatrixBuilder::LeftRightColIndex(colCount, left))
      );
      MATHICGB_ASSERT(inserted.second);
      MATHICGB_ASSERT(inserted.first.first != 0);
      auto p(inserted.first);
      toMonomial.back() = copied;

      MATHICGB_ASSERT(ring.monomialEqualHintTrue(copied, p.second));
      MATHICGB_ASSERT(*p.first == QuadMatrixBuilder::LeftRightColIndex(colCount, left));
      return ::std::make_pair(*p.first, p.second);
    } catch (...) {
      toMonomial.pop_back();
      ring.freeMonomial(copied);
      throw;
    }
  }
}

::std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial>
QuadMatrixBuilder::createColumnLeft(
  const_monomial monomialToBeCopied
) {
  return createCol
    (monomialToBeCopied,
     mTopLeft,
     mBottomLeft,
     mMonomialsLeft,
     mMonomialToCol,
     ring(),
     true);
}

::std::pair<QuadMatrixBuilder::LeftRightColIndex, ConstMonomial>
QuadMatrixBuilder::createColumnRight(
  const_monomial monomialToBeCopied
) {
  return createCol
    (monomialToBeCopied,
     mTopRight,
     mBottomRight,
     mMonomialsRight,
     mMonomialToCol,
     ring(),
     false);
}

QuadMatrix QuadMatrixBuilder::buildMatrixAndClear() {
  QuadMatrix out;

  mTopLeft.swap(out.topLeft);
  mTopRight.swap(out.topRight);
  mBottomLeft.swap(out.bottomLeft);
  mBottomRight.swap(out.bottomRight);
  out.ring = &ring();

  mTopLeft.clear();
  mTopRight.clear();
  mBottomLeft.clear();
  mBottomRight.clear();

  MATHICGB_ASSERT(out.debugAssertValid());
  return ::std::move(out);
}

MATHICGB_NAMESPACE_END
