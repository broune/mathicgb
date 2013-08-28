// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "DivisorLookup.hpp"

#include "SigPolyBasis.hpp"
#include "DivLookup.hpp"
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

namespace {
  template<
    template<class> class BaseLookup,
    bool AllowRemovals,
    bool UseDivMask
  >
  std::unique_ptr<DivisorLookup> create(
    const DivisorLookup::Monoid& monoid,
    int type,
    bool preferSparseReducers
  ) {
    auto p = new DivLookup<BaseLookup, AllowRemovals, UseDivMask>
      (monoid, type, preferSparseReducers);
    return std::unique_ptr<DivisorLookup>(p);
  }

  std::unique_ptr<DivisorLookup> createGeneral(
    const DivisorLookup::Monoid& monoid,
    int type,
    bool preferSparseReducers,
    bool allowRemovals
  ) {
    const bool sparse = preferSparseReducers;
    switch (type) {
    case 1:
      if (allowRemovals)
        return create<mathic::DivList, true, true>(monoid, type, sparse);
      else
        return create<mathic::DivList, false, true>(monoid, type, sparse);

    case 2:
      if (allowRemovals)
        return create<mathic::KDTree, true, true>(monoid, type, sparse);
      else
        return create<mathic::KDTree, false, true>(monoid, type, sparse);

    case 3:
      if (allowRemovals)
        return create<mathic::DivList, true, false>(monoid, type, sparse);
      else
        return create<mathic::DivList, false, false>(monoid, type, sparse);

    case 4:
      if (allowRemovals)
        return create<mathic::KDTree, true, false>(monoid, type, sparse);
      else
        return create<mathic::KDTree, false, false>(monoid, type, sparse);

    default:
      MATHICGB_ASSERT_NO_ASSUME(false);
      throw std::runtime_error("Unknown code for monomial data structure");
    }
  }

  class ConcreteFactory : public DivisorLookup::Factory {
  public:
    ConcreteFactory(const Monoid& monoid, int type): 
      mMonoid(monoid),
      mType(type)
    {}

    virtual std::unique_ptr<DivisorLookup> create(
      bool preferSparseReducers,
      bool allowRemovals
    ) const {
      return createGeneral
        (mMonoid, mType, preferSparseReducers, allowRemovals);
    }

  private:
    const Monoid& mMonoid;
    const int mType;
  };
}

std::unique_ptr<DivisorLookup::Factory> DivisorLookup::makeFactory(
  const PolyRing& ring,
  int type
) {
  return std::unique_ptr<Factory>(new ConcreteFactory(ring.monoid(), type));
}

void DivisorLookup::displayDivisorLookupTypes(std::ostream &o)
{
  o << "Divisor Lookup Types:" << std::endl;
  o << "  1   divlist+divmask" << std::endl;
  o << "  2   kdtree+divmask" << std::endl;
  o << "  3   divlist" << std::endl;
  o << "  4   kdtree" << std::endl;
}

MATHICGB_NAMESPACE_END
