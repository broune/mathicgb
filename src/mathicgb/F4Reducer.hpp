#ifndef MATHICGB_F4_REDUCER_GUARD
#define MATHICGB_F4_REDUCER_GUARD

#include "Reducer.hpp"
#include "PolyRing.hpp"
#include <string>
class QuadMatrix;

class F4Reducer : public Reducer {
public:
  enum Type {
    OldType,
    NewType
  };

  F4Reducer(const PolyRing& ring, Type type);

  virtual size_t preferredSetSize() const;

  /// Store all future matrices to file-1.mat, file-2.mat and so on.
  /// Matrices with less than minEntries non-zero entries are not stored.
  /// If file is an empty string then no matrices are stored. If this method
  /// is never called then no matrices are stored.
  void writeMatricesTo(std::string file, size_t minEntries);

  virtual std::unique_ptr<Poly> classicReduce
    (const Poly& poly, const PolyBasis& basis);

  virtual std::unique_ptr<Poly> classicTailReduce
    (const Poly& poly, const PolyBasis& basis);

  virtual std::unique_ptr<Poly> classicReduceSPoly
    (const Poly& a, const Poly& b, const PolyBasis& basis);

  virtual void classicReduceSPolySet(
    std::vector<std::pair<size_t, size_t> >& spairs,
    const PolyBasis& basis,
    std::vector<std::unique_ptr<Poly> >& reducedOut
  );

  virtual void classicReducePolySet(
    const std::vector<std::unique_ptr<Poly> >& polys,
    const PolyBasis& basis,
    std::vector<std::unique_ptr<Poly> >& reducedOut
  );

  virtual Poly* regularReduce(
    const_monomial sig,
    const_monomial multiple,
    size_t basisElement,
    const SigPolyBasis& basis
  );

  virtual void setMemoryQuantum(size_t quantum);

  virtual std::string description() const;
  virtual size_t getMemoryUse() const;

private:
  void saveMatrix(const QuadMatrix& matrix);

  Type mType;
  std::unique_ptr<Reducer> mFallback;
  const PolyRing& mRing;
  size_t mMemoryQuantum;
  std::string mStoreToFile; /// stem of file names to save matrices to
  size_t mMinEntryCountForStore; /// don't save matrices with fewer entries
  size_t mMatrixSaveCount; // how many matrices have been saved
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
