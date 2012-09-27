#ifndef F4_MATRIX_BUILDER_GUARD
#define F4_MATRIX_BUILDER_GUARD

#include "QuadMatrixBuilder.hpp"
#include "Poly.hpp"
#include "PolyRing.hpp"
#include "PolyBasis.hpp"
#include "QuadMatrix.hpp"
#include <vector>

/** Class for constructing an F4 matrix. This class is reponsible for
  figuring out what matrix to build and then it uses QuadMatrixBuilder
  to create that matrix.

  @todo: this class does not offer exception guarantees. It's just not
  very workable without an RAII monomial handle, so add one of those
  before fixing this. */
class F4MatrixBuilder {
private:
  typedef QuadMatrixBuilder::ColIndex ColIndex;
  typedef QuadMatrixBuilder::LeftRightColIndex LeftRightColIndex;

public:
  F4MatrixBuilder(const PolyBasis& basis);

  /** Schedules two rows to be added to the matrix whose linear span
    includes the S-polynomial between polyA and polyB. More precisely,
    the two rows represent (B:A)*polyA and (A:B)*polyB where
    A=lead(polyA) and B=lead(polyB). */
  void addTwoRowsForSPairToMatrix(const Poly& polyA, const Poly& polyB);

  /** Schedules a row representing multiple*poly to be added to the
      matrix. No ownership is taken, but poly must remain valid until
      the matrix is constructed. multiple is copied so there is no
      requirement there. */
  void addRowToMatrix(const_monomial multiple, const Poly& poly);

  /** Builds an F4 matrix to the specifications given. Also clears the
    information in this object.

    The right columns are ordered by decreasing monomial of that
    column according to the order from the basis. The left columns are
    order in some way so that the first entry in each row (the pivot)
    has a lower index than any other entries in that row.

    The matrix contains a reducer/pivot for every monomial that can be
    reduced by the basis and that is present in the matrix. Note that
    the client-added rows also count as reducers so their lead terms
    will not get another reducer added automatically -- specifically,
    adding an S-polynomial will not do what you want because its lead
    term may have no reducer in the matrix other than itself. Instead,
    add the two polynomials that you would have subtracted from each
    other to form the S-polynomial - or even better call the method
    that adds an S-pair for you. */
  void buildMatrixAndClear(QuadMatrix& matrix);

  const PolyRing& ring() const {return mBuilder.ring();}

private:
  /** Returns a left or right column that represents mono. Creates a
      new column and schedules a new row to reduce that column if necessary. */
  LeftRightColIndex createOrFindColumnOf(const_monomial mono);

  void appendRowTop(const_monomial multiple, const Poly& poly);
  void appendRowBottom(const_monomial multiple, const Poly& poly);

  /// Represents the task of adding a row representing poly*multiple.
  struct RowTask {
    bool useAsReducer; // if true: put in top part of matrix
    const Poly* poly;
    monomial multiple;
  };
  std::vector<RowTask> mTodo;
  const PolyBasis& mBasis;
  QuadMatrixBuilder mBuilder;
};

#endif
