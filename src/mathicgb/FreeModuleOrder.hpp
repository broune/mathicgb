// Copyright 2011 Michael E. Stillman

#ifndef _free_module_order_h_
#define _free_module_order_h_

#include "PolyRing.hpp"
#include "SigSPairQueue.hpp"

class PolyRing;
class GroebnerBasis;
class SigSPairQueue;
typedef int FreeModuleOrderType;

class FreeModuleOrder
{
public:
  FreeModuleOrder() {}
  virtual ~FreeModuleOrder() {}

  // returns LT, EQ, or GT, depending on sig ? sig2.
  virtual int signatureCompare(const_monomial sig, const_monomial sig2) const = 0;

  // compares sig vs (m2*sig)
  virtual int signatureCompare(
    const_monomial sig,
    const_monomial m2,
    const_monomial sig2
  ) const = 0;

  // Sorts in ascending order of signature.
  virtual void sortSignatures(std::vector<PreSPair>& pairs) const = 0;

  virtual void getStats(size_t& comparisons, size_t& preComparisons) const = 0;

  virtual std::string description() const = 0;

  virtual std::unique_ptr<SigSPairQueue>
  createSigSPairQueue(GroebnerBasis const& basis) const = 0;

  /// @todo: We need at least an enum to make this clearer
  static std::unique_ptr<FreeModuleOrder> makeOrder(FreeModuleOrderType type, const PolyRing& ring);

  static void displayOrderTypes(std::ostream &o);

};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
