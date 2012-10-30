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
  typedef QuadMatrixBuilder::Scalar Scalar;

public:
  F4MatrixBuilder(const PolyBasis& basis, int threadCount);

  /** Schedules a row representing the S-polynomial between polyA and
    polyB to be added to the matrix. No ownership is taken, but polyA
    and polyB must remain valid until the matrix is constructed.

    Currently, the two monomials must be monic, though this is just
    because they always happen to be monic so there was no reason to
    support the non-monic case. */
  void addSPolynomialToMatrix(const Poly& polyA, const Poly& polyB);

  /** Schedules a row representing multiple*poly to be added to the
    matrix. No ownership is taken, but poly must remain valid until
    the matrix is constructed. multiple is copied so it need not
    remain valid. */
  void addPolynomialToMatrix(const_monomial multiple, const Poly& poly);

  /// as the overload with a multiple, just letting multiple be the
  /// identity.
  void addPolynomialToMatrix(const Poly& poly);

  /** Builds an F4 matrix to the specifications given. Also clears the
    information in this object.

    The right columns are ordered by decreasing monomial of that
    column according to the order from the basis. The left columns are
    ordered in some way so that the first entry in each top row (the
    pivot) has a lower index than any other entries in that row.

    The matrix contains a reducer/pivot for every monomial that can be
    reduced by the basis and that is present in the matrix. There is
    no guarantee that the bottom part of the matrix contains rows that
    exactly correspond to the polynomials that have been scheduled to
    be added to the matrix. It is only guaranteed that the matrix has
    the same row-space as though that had been the case. */
  void buildMatrixAndClear(QuadMatrix& matrix);

  const PolyRing& ring() const {return mBuilder.ring();}

private:
  /** Creates a column with monomial label x and schedules a new row to
    reduce that column if possible. Here x is monoA if monoB is
    null and otherwise x is the product of monoA and monoB. */
  MATHICGB_NO_INLINE LeftRightColIndex createColumn
    (QuadMatrixBuilder& builder, const_monomial monoA, const_monomial monoB);

  /// Represents the task of adding a row to the matrix. If sPairPoly is null
  /// then the row to add is multiply * poly. Otherwise, the row to add is
  ///   multiply * poly - sPairMultiply * sPairPoly
  /// where sPairMultiply makes the lead terms cancel.
  struct RowTask {
    bool addToTop; // add the row to the bottom if false
    const Poly* poly;
    monomial multiply;
    const Poly* sPairPoly;
    monomial sPairMultiply;
  };

  void appendRowTop(
    const_monomial multiple,
    const Poly& poly,
    QuadMatrixBuilder& builder);
  void appendRowBottom(const RowTask& task, QuadMatrixBuilder& builder);
  void appendRowBottom(
    const_monomial multiple,
    bool negate,
    Poly::const_iterator begin,
    Poly::const_iterator end,
    QuadMatrixBuilder& builder
  );

  MATHICGB_INLINE LeftRightColIndex findOrCreateColumn
    (const_monomial monoA, const_monomial monoB, QuadMatrixBuilder& builder);


  const int mThreadCount;
  monomial mTmp;
  std::vector<RowTask> mTodo;
  const PolyBasis& mBasis;
  QuadMatrixBuilder mBuilder;
};

#endif
