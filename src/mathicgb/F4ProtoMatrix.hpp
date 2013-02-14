#ifndef MATHICGB_F4_PROTO_MATRIX_GUARD
#define MATHICGB_F4_PROTO_MATRIX_GUARD

#include "PolyRing.hpp"
#include "SparseMatrix.hpp"
#include "Poly.hpp"

class F4ProtoMatrix {
public:
  typedef uint32 RowIndex;
  typedef uint32 ColIndex;
  typedef coefficient ExternalScalar;
  typedef SparseMatrix::Scalar Scalar;

  struct Row {
    const ColIndex* indices;
    const Scalar* scalars;
    const ExternalScalar* externalScalars;
    ColIndex entryCount;
  };

  RowIndex rowCount() const {return static_cast<RowIndex>(mRows.size());}

  Row row(const RowIndex row) const {
    MATHICGB_ASSERT(row < mRows.size());
    const auto& r = mRows[row];
    Row rr;
    rr.indices = mIndices.data() + r.indicesBegin;
    rr.entryCount = r.entryCount;
    if (r.externalScalars == 0) {
      rr.scalars = mScalars.data() + r.scalarsBegin;
      rr.externalScalars = 0;
    } else {
      rr.scalars = 0;
      rr.externalScalars = r.externalScalars;
    }
    return rr;
  }

  ColIndex* makeRowWithTheseScalars(const Poly& scalars);

  std::pair<ColIndex*, Scalar*> makeRow(ColIndex entryCount);

  void removeLastEntries(const RowIndex row, const ColIndex count);

private:
  struct InternalRow {
    size_t indicesBegin;
    size_t scalarsBegin;
    ColIndex entryCount;
    const ExternalScalar* externalScalars;
  };

  std::vector<ColIndex> mIndices;
  std::vector<Scalar> mScalars;
  std::vector<InternalRow> mRows;
};

#endif
