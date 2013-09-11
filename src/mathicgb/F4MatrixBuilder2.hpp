// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_F4_MATRIX_BUILDER_2_GUARD
#define MATHICGB_F4_MATRIX_BUILDER_2_GUARD

#include "SparseMatrix.hpp"
#include "Poly.hpp"
#include "PolyRing.hpp"
#include "PolyBasis.hpp"
#include "QuadMatrix.hpp"
#include "MonomialMap.hpp"
#include "F4ProtoMatrix.hpp"
#include "mtbb.hpp"
#include <vector>

MATHICGB_NAMESPACE_BEGIN

/// Class for constructing an F4 matrix.
///
/// @todo: this class does not offer exception guarantees. It's just not
/// very workable without an RAII monomial handle or a scope exit
/// functionality, so add one of those before fixing this.
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

  /// Schedules a row representing the S-polynomial between polyA and
  /// polyB to be added to the matrix. No ownership is taken, but polyA
  /// and polyB must remain valid until the matrix is constructed.
  ///
  /// Currently, the two monomials must be monic, though this is just
  /// because they happen always to be monic so there was no reason to
  /// support the non-monic case.
  void addSPolynomialToMatrix(const Poly& polyA, const Poly& polyB);

  /// Schedules a row representing multiple*poly to be added to the
  /// matrix. No ownership is taken, but poly must remain valid until
  /// the matrix is constructed. multiple is copied, so it need not
  /// remain valid.
  void addPolynomialToMatrix(const_monomial multiple, const Poly& poly);

  /// As the overload with a multiple, where the multiple is 1.
  void addPolynomialToMatrix(const Poly& poly);

  /// Builds an F4 matrix to the specifications given. Also clears the
  /// information in this object.
  ///
  /// The right columns are in order of strictly decreasing monomial.
  /// The left columns are ordered in some way so that the leading non-zero
  /// entry in each top row has the maximal column monomial out of all
  /// non-zero entries in that row.
  ///
  /// The monomials that can be reduced by some element of the basis go on
  /// the left while the remaining monomials go on the right. The upper left
  /// matrix is upper triangular, thus having a reducer/pivot row for every
  /// column.
  ///
  /// There is no guarantee that the bottom part of the matrix contains rows
  /// that exactly correspond to the polynomials that have been scheduled to
  /// be added to the matrix. It is only guaranteed that the whole matrix has
  /// the same row-space as though that had been the case.
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
    monomial desiredLead; // multiply a monomial onto poly to get this lead
    const Poly* poly;
    const Poly* sPairPoly;
  };
  typedef mgb::mtbb::parallel_do_feeder<RowTask> TaskFeeder;

  /// Creates a column with monomial label monoA * monoB and schedules a new
  /// row to reduce that column if possible. If such a column already
  /// exists, then a new column is not inserted. In either case, returns
  /// the column index and column monomial corresponding to monoA * monoB.
  ///
  /// createColumn can be used simply to search for an existing column, but
  /// since createColumn incurs locking overhead, this is not a good idea.
  /// Note that createColumn has to work correctly for pre-existing columns
  /// because the only way to be *certain* that no other thread has inserted
  /// the column of interest is to grab a lock, and the lock being grabbed
  /// is being grabbed inside createColumn.
  MATHICGB_NO_INLINE
  std::pair<ColIndex, ConstMonomial> createColumn(
    const_monomial monoA,
    const_monomial monoB,
    TaskFeeder& feeder
  );

  /// Append multiple * poly to block, creating new columns as necessary.
  void appendRow(
    const_monomial multiple,
    const Poly& poly,
    F4ProtoMatrix& block,
    TaskFeeder& feeder
  );

  /// Append poly*multiply - sPairPoly*sPairMultiply to block, creating new
  /// columns as necessary.
  void appendRowSPair(
    const Poly* poly,
    monomial multiply,
    const Poly* sPairPoly,
    monomial sPairMultiply,
    F4ProtoMatrix& block,
    TaskFeeder& feeder
  );

  /// As createColumn, except with much better performance in the common
  /// case that the column for monoA * monoB already exists. In particular,
  /// no lock is grabbed in that case.
  MATHICGB_NO_INLINE
  std::pair<ColIndex, ConstMonomial> findOrCreateColumn(
    const_monomial monoA,
    const_monomial monoB,
    TaskFeeder& feeder
  );

  /// As the overload that does not take a ColReader parameter, except with
  /// better performance in the common case that the column already exists
  /// and colMap is up-to-date.
  MATHICGB_INLINE
  std::pair<ColIndex, ConstMonomial> findOrCreateColumn(
    const_monomial monoA,
    const_monomial monoB,
    const ColReader& colMap,
    TaskFeeder& feeder
  );

  /// The split into left and right columns is not done until the whole matrix
  /// has been constructed. This vector keeps track of which side each column
  /// should go to once we do the split. char is used in place of bool because
  /// the specialized bool would just be slower for this use case. See
  /// http://isocpp.org/blog/2012/11/on-vectorbool .
  std::vector<char> mIsColumnToLeft;

  /// How much memory to allocate every time more memory is needed.
  const size_t mMemoryQuantum;

  /// If you want to modify the columns, you need to grab this lock first.
  mgb::mtbb::mutex mCreateColumnLock;

  /// A monomial for temporary scratch calculations. Protected by
  /// mCreateColumnLock.
  monomial mTmp;

  /// The basis that supplies reducers.
  const PolyBasis& mBasis;

  /// Mapping from monomials to column indices.
  Map mMap;

  /// Stores the rows that have been scheduled to be added.
  std::vector<RowTask> mTodo;
};

MATHICGB_NAMESPACE_END
#endif
