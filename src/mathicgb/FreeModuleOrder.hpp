// Copyright 2011 Michael E. Stillman

#ifndef _free_module_order_h_
#define _free_module_order_h_

#include "PolyRing.hpp"
#include "PairTriangle.hpp"

class PolyRing;
class SPairQueue;
class Ideal;
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

  // Sorts in ascending order of signature. May alter the signatures,
  // so do not use them after calling this method other than to free them.
  virtual void destructiveSort(std::vector<PreSPair>& pairs) const = 0;

  // You must use this method to inform the order when a
  // new basis element has been added.
  virtual void appendBasisElement(const_monomial m) = 0;

  virtual void getStats(size_t& comparisons, size_t& preComparisons) const = 0;

  virtual std::string description() const = 0;

  virtual std::auto_ptr<SPairQueue> makeQueue(size_t type) const = 0;

  static FreeModuleOrder* makeOrder(FreeModuleOrderType type, const Ideal* I);

  static void displayOrderTypes(std::ostream &o);
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
