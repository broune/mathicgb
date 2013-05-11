// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include "FreeModuleOrder.hpp"

#include "Poly.hpp"
#include "SigSPairQueue.hpp"
#include "GroebnerBasis.hpp"
#include "PolyRing.hpp"
#include <mathic.h>
#include <iostream>
#include <algorithm>
#include <limits>
#include <stdexcept>

class ConcreteOrder : public FreeModuleOrder {
public:
  ConcreteOrder(const PolyRing& ring): mRing(ring) {}

  virtual int signatureCompare(const_monomial sigA, const_monomial sigB) const {
    return mRing.monomialCompare(sigA, sigB);
  }

  virtual int signatureCompare(
    const_monomial sigA,
    const_monomial monoB,
    const_monomial sigB
  ) const {
    return mRing.monomialCompare(sigA, monoB, sigB);
  }

  virtual std::string description() const {
    return "todo";
  }

private:
  const PolyRing& mRing;
};

std::unique_ptr<FreeModuleOrder> FreeModuleOrder::makeOrder(FreeModuleOrderType type, const PolyRing& ring)
{
  return make_unique<ConcreteOrder>(ring);
}
