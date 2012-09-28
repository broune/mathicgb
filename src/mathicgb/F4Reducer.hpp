#ifndef MATHICGB_F4_REDUCER_GUARD
#define MATHICGB_F4_REDUCER_GUARD

#include "Reducer.hpp"
#include "PolyRing.hpp"

class F4Reducer : public Reducer {
public:
  F4Reducer(const PolyRing& ring,
            std::unique_ptr<Reducer> fallback);

  virtual std::unique_ptr<Poly> classicReduce
  (const Poly& poly, const PolyBasis& basis);

  virtual std::unique_ptr<Poly> classicTailReduce
  (const Poly& poly, const PolyBasis& basis);

  virtual std::unique_ptr<Poly> classicReduceSPoly
  (const Poly& a, const Poly& b, const PolyBasis& basis);

  virtual void classicReduceSPolyGroup
  (std::vector<std::pair<size_t, size_t> >& spairs,
   const PolyBasis& basis,
   std::vector<std::unique_ptr<Poly> >& reducedOut);

  virtual Poly* regularReduce(
    const_monomial sig,
    const_monomial multiple,
    size_t basisElement,
    const GroebnerBasis& basis);

  virtual void setThreadCount(size_t threadCount);

  virtual std::string description() const;
  virtual size_t getMemoryUse() const;

private:
  std::unique_ptr<Reducer> mFallback;
  const PolyRing& mRing;
  size_t mThreadCount;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
