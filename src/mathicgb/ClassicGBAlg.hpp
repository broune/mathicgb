// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_CLASSIC_GB_ALG_GUARD
#define MATHICGB_CLASSIC_GB_ALG_GUARD

#include "Reducer.hpp"
#include "SPairs.hpp"
#include "PolyBasis.hpp"
#include <mathic.h>
#include <memory>
#include <ostream>
#include <vector>

MATHICGB_NAMESPACE_BEGIN

class Basis;

/// Calculates a classic Grobner basis using Buchberger's algorithm.
class ClassicGBAlg {
public:
  ClassicGBAlg(
    const Basis& basis,
    Reducer& reducer,
    int divisorLookupType,
    bool preferSparseReducers,
    size_t queueType
  );

  // Replaces the current basis with a Grobner basis of the same ideal.
  void computeGrobnerBasis();

  // How many S-pairs were not eliminated before reduction of the
  // corresponding S-polynomial.
  unsigned long long sPolyReductionCount() const {return mSPolyReductionCount;}

  // Returns the current basis.
  PolyBasis& basis() {return mBasis;}

  // Shows statistics on what the algorithm has done.
  void printStats(std::ostream& out) const;

  void printMemoryUse(std::ostream& out) const;

  size_t getMemoryUse() const;

  void setBreakAfter(unsigned int elements) {
    mBreakAfter = elements;
  }

  void setPrintInterval(unsigned int reductions) {
    mPrintInterval = reductions;
  }

  /// A value of zero means to let the algorithm decide a reasonable
  /// value based on the other settings.
  void setSPairGroupSize(unsigned int groupSize);

  void setReducerMemoryQuantum(size_t memoryQuantum) {
    mReducer.setMemoryQuantum(memoryQuantum);
  }

  void setUseAutoTopReduction(bool value) {
    mUseAutoTopReduction = value;
  }

  void setUseAutoTailReduction(bool value) {
    mUseAutoTailReduction = value;
  }

  class Callback {
  public:
    /// Stop the computation if call return false.
    virtual bool call() = 0;
  };
  /// callback is called every once in a while and then it has the
  /// option of stopping the computation. callback can be null, in
  /// which case no call is made and the computation continues.
  void setCallback(Callback* callback) {mCallback = callback;}

private:
  Callback* mCallback;
  unsigned int mBreakAfter;
  unsigned int mPrintInterval;
  unsigned int mSPairGroupSize;
  bool mUseAutoTopReduction;
  bool mUseAutoTailReduction;

  // Perform a step of the algorithm.
  void step();

  void autoTailReduce();

  void insertReducedPoly(std::unique_ptr<Poly> poly);

  // clears polynomials.
  void insertPolys(std::vector<std::unique_ptr<Poly> >& polynomials);

  const PolyRing& mRing;
  Reducer& mReducer;
  PolyBasis mBasis;
  SPairs mSPairs;
  mic::Timer mTimer;
  unsigned long long mSPolyReductionCount;
};

MATHICGB_NAMESPACE_END
#endif
