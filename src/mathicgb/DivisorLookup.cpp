// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "DivisorLookup.hpp"

#include "SigPolyBasis.hpp"
#include "DivLookup.hpp"
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

namespace {
  struct DefaultParams {
    static bool const minimizeOnInsert = false;
    static bool const sortOnInsert = false;
    static bool const useDivisorCache = true;
    static double const rebuildRatio;
    static size_t const minRebuildRatio = 50;
  };
  double const DefaultParams::rebuildRatio = 0.5;

  class DivListFactory : public DivisorLookup::Factory {
  public:
    DivListFactory(const PolyRing& ring, bool useDivMask): mRing(ring), mUseDivMask(useDivMask) {}
    virtual std::unique_ptr<DivisorLookup> create(
      bool preferSparseReducers,
      bool allowRemovals
    ) const {
      if (mUseDivMask)
        return createIt<true>(preferSparseReducers);
      else
        return createIt<false>(preferSparseReducers);
    }

  private:
    template<bool UseDivMask>
    std::unique_ptr<DivisorLookup> createIt(bool preferSparseReducers) const {
      typedef DivLookupConfiguration<true, UseDivMask> Configuration;
      bool useDM = UseDivMask;
      Configuration configuration(
        mRing,
        DefaultParams::minimizeOnInsert,
        DefaultParams::sortOnInsert,
        DefaultParams::useDivisorCache,
        DefaultParams::rebuildRatio,
        DefaultParams::minRebuildRatio,
        (useDM ? 1 : 3),
        preferSparseReducers);
      return std::unique_ptr<DivisorLookup>
        (new DivLookup<mathic::DivList<Configuration> >(configuration));
    }

    const PolyRing& mRing;
    bool mUseDivMask;
  };

  class KDTreeFactory : public DivisorLookup::Factory {
  public:
    KDTreeFactory(const PolyRing& ring, bool useDivMask): mRing(ring), mUseDivMask(useDivMask) {}
    virtual std::unique_ptr<DivisorLookup> create(
      bool preferSparseReducers,
      bool allowRemovals
    ) const {
      if (allowRemovals) {
        if (mUseDivMask)
          return createAllowRemovals<true,true>(preferSparseReducers);
        else
          return createAllowRemovals<true,false>(preferSparseReducers);
      } else {
        if (mUseDivMask)
          return createAllowRemovals<false,true>(preferSparseReducers);
        else
          return createAllowRemovals<false,false>(preferSparseReducers);
      }
    }

  private:
    template<bool AllowRemovals, bool UseDivMask>
    std::unique_ptr<DivisorLookup> createAllowRemovals(
      bool preferSparseReducers
    ) const {
      typedef DivLookupConfiguration<AllowRemovals, UseDivMask> Configuration;
      bool useDM = UseDivMask;
      Configuration configuration(
        mRing,
        DefaultParams::minimizeOnInsert,
        DefaultParams::sortOnInsert,
        DefaultParams::useDivisorCache,
        DefaultParams::rebuildRatio,
        DefaultParams::minRebuildRatio,
        (useDM ? 2 : 4),
        preferSparseReducers);
      return std::unique_ptr<DivisorLookup>
        (new DivLookup<mathic::KDTree<Configuration> >(configuration));
    }
    const PolyRing& mRing;
    bool mUseDivMask;
  };
}

std::unique_ptr<DivisorLookup::Factory> DivisorLookup::makeFactory(
  const PolyRing& ring,
  int type
) {
  if (type == 1)
    return std::unique_ptr<Factory>(new DivListFactory(ring, true));
  else if (type == 2)
    return std::unique_ptr<Factory>(new KDTreeFactory(ring, true));
  if (type == 3)
    return std::unique_ptr<Factory>(new DivListFactory(ring, false));
  else if (type == 4)
    return std::unique_ptr<Factory>(new KDTreeFactory(ring, false));
  else if (type == 0)
    throw std::runtime_error("Divisor lookup 0 (DivisorLookupGB) disabled.");
  else
    throw std::runtime_error("Unknown divisor lookup code.");
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
