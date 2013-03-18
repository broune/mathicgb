#ifndef F4_MATRIX_BUILDER_2_GUARD
#define F4_MATRIX_BUILDER_2_GUARD

#include "SparseMatrix.hpp"
#include "Poly.hpp"
#include "PolyRing.hpp"
#include "PolyBasis.hpp"
#include "QuadMatrix.hpp"
#include "MonomialMap.hpp"
#include "F4ProtoMatrix.hpp"
#include "mtbb.hpp"
#include <vector>

/** Class for constructing an F4 matrix.

  @todo: this class does not offer exception guarantees. It's just not
  very workable without an RAII monomial handle or a scope exit functionality,
  so add one of those before fixing this. */
class F4MatrixBuilder2 {
private:
  typedef SparseMatrix::ColIndex ColIndex;
  typedef SparseMatrix::Scalar Scalar;
  typedef MonomialMap<ColIndex> Map;
  typedef SparseMatrix::RowIndex RowIndex;

public:
  /// memoryQuantum is how much to increase the memory size by each time the
  /// current amount of memory is exhausted. A value of 0 indicates to start
  /// small and double the quantum at each exhaustion.
  F4MatrixBuilder2(const PolyBasis& basis, size_t memoryQuantum = 0);

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

    The right columns are ordered by decreasing monomial of each
    column according to the order from the basis. The left columns are
    ordered in some way so that the first entry in each top row (the
    pivot) has a lower index than any other entries in that row.

    The matrix contains a reducer/pivot for every monomial that can be
    reduced by the basis and that is present in the matrix. There is
    no guarantee that the bottom part of the matrix contains rows that
    exactly correspond to the polynomials that have been scheduled to
    be added to the matrix. It is only guaranteed that the whole matrix has
    the same row-space as though that had been the case. */
  void buildMatrixAndClear(QuadMatrix& matrix);

  const PolyRing& ring() const {return mBasis.ring();}

private:
  typedef const Map::Reader ColReader;
  typedef std::vector<monomial> Monomials;

  /// Represents the task of adding a row to the matrix. If sPairPoly is null
  /// then the row to add is multiply * poly. Otherwise, the row to add is
  ///   multiply * poly - sPairMultiply * sPairPoly
  /// where multiply and sPairMultiply are such that the leading terms become
  /// desiredLead.
  struct RowTask {
    monomial desiredLead; // multiply monomial onto poly to get this lead
    const Poly* poly;
    const Poly* sPairPoly;
  };
  typedef mgb::tbb::parallel_do_feeder<RowTask> TaskFeeder;

  /// Creates a column with monomial label x and schedules a new row to
  /// reduce that column if possible. Here x is monoA if monoB is
  /// null and otherwise x is the product of monoA and monoB.
  MATHICGB_NO_INLINE
  std::pair<ColIndex, ConstMonomial> createColumn(
    const_monomial monoA,
    const_monomial monoB,
    TaskFeeder& feeder
  );

  void appendRow(
    const_monomial multiple,
    const Poly& poly,
    F4ProtoMatrix& block,
    TaskFeeder& feeder
  );
  void appendRowSPair(
    const Poly* poly,
    monomial multiply,
    const Poly* sPairPoly,
    monomial sPairMultiply,
    F4ProtoMatrix& block,
    TaskFeeder& feeder
  );

  MATHICGB_NO_INLINE
  std::pair<ColIndex, ConstMonomial> findOrCreateColumn(
    const_monomial monoA,
    const_monomial monoB,
    TaskFeeder& feeder
  );
  
  MATHICGB_INLINE
  std::pair<ColIndex, ConstMonomial> findOrCreateColumn(
    const_monomial monoA,
    const_monomial monoB,
    const ColReader& colMap,
    TaskFeeder& feeder
  );

  MATHICGB_NO_INLINE
  void createTwoColumns(
    const const_monomial monoA1,
    const const_monomial monoA2,
    const const_monomial monoB,
    TaskFeeder& feeder
  );

  std::vector<char> mIsColumnToLeft;
  const size_t mMemoryQuantum;
  mgb::tbb::mutex mCreateColumnLock;
  monomial mTmp;
  const PolyBasis& mBasis;
  Map mMap;
  std::vector<RowTask> mTodo;
};

#endif
