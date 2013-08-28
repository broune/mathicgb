// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "MonoLookup.hpp"

#include "SigPolyBasis.hpp"
#include "StaticMonoLookup.hpp"
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

namespace {
  template<
    template<class> class BaseLookupTemplate,
    bool AllowRemovals,
    bool UseDivMask
  >
  class ConcreteMonoLookup : public MonoLookup {
  public:
    ConcreteMonoLookup(
      const Monoid& monoid,
      int type,
      bool preferSparseReducers
    ):
      mLookup(monoid),
      mType(type),
      mPreferSparseReducers(preferSparseReducers),
      mBasis(0),
      mSigBasis(0)
    {}

    const Monoid& monoid() const {return mLookup.monoid();}
    bool preferSparseReducers() const {return mPreferSparseReducers;}

    // *** Virtual interface follows

    virtual void setBasis(const PolyBasis& basis) {
      if (mBasis == &basis)
        return;
      MATHICGB_ASSERT(mBasis == 0);
      MATHICGB_ASSERT(monoid() == basis.ring().monoid());
      mBasis = &basis;
    }

    virtual void setSigBasis(const SigPolyBasis& sigBasis) {
      if (mSigBasis == &sigBasis)
        return;
      MATHICGB_ASSERT(mSigBasis == 0);
      MATHICGB_ASSERT(mBasis == 0 || mBasis == &sigBasis.basis());
      MATHICGB_ASSERT(monoid() == sigBasis.basis().ring().monoid());
      mSigBasis = &sigBasis;
      setBasis(sigBasis.basis());
    }

    const SigPolyBasis& sigBasis() const {
      MATHICGB_ASSERT(mSigBasis != 0);
      return *mSigBasis;
    }

    const PolyBasis& basis() const {
      MATHICGB_ASSERT(mBasis != 0);
      return *mBasis;
    }

    virtual void insert(ConstMonoRef mono, size_t index) {
      mLookup.insert(mono, index);
    }

    virtual size_t regularReducer(ConstMonoRef sig, ConstMonoRef mono) const {
      return mLookup.regularReducer
        (sig, mono, sigBasis(), preferSparseReducers());
    }

    virtual size_t classicReducer(ConstMonoRef mono) const {
      return mLookup.classicReducer(mono, basis(), preferSparseReducers());
    }

    virtual std::string getName() const {return mLookup.getName();}

    virtual size_t getMemoryUse() const {return mLookup.getMemoryUse();}

    virtual size_t highBaseDivisor(size_t newGenerator) const {
      return mLookup.highBaseDivisor(newGenerator, sigBasis());
    }
      
    virtual void lowBaseDivisors(
      std::vector<size_t>& divisors,
      size_t maxDivisors,
      size_t newGenerator
    ) const {
      return mLookup.lowBaseDivisors
        (divisors, maxDivisors, newGenerator, sigBasis());
    }
    virtual size_t minimalLeadInSig(ConstMonoRef sig) const {
      return mLookup.minimalLeadInSig(sig, sigBasis());
    }

    virtual int type() const {return mType;}

    virtual void multiples(ConstMonoRef mono, EntryOutput& consumer) const {
      mLookup.multiples(mono, consumer);
    }

    virtual size_t divisor(ConstMonoRef mono) const {
      return mLookup.divisor(mono);
    }

    virtual void divisors(ConstMonoRef mono, EntryOutput& consumer) const {
      mLookup.divisors(mono, consumer);
    }

    virtual void removeMultiples(ConstMonoRef mono) {
      mLookup.removeMultiples(mono);
    }

    virtual void remove(ConstMonoRef mono) {
      return mLookup.remove(mono);
    }

    virtual size_t size() const {return mLookup.size();}

  private:
    StaticMonoLookup<BaseLookupTemplate, AllowRemovals, UseDivMask> mLookup;
    const int mType;
    const bool mPreferSparseReducers;
    PolyBasis const* mBasis;
    SigPolyBasis const* mSigBasis;
  };

  template<
    template<class> class BaseLookup,
    bool AllowRemovals,
    bool UseDivMask
  >
  std::unique_ptr<MonoLookup> create(
    const MonoLookup::Monoid& monoid,
    int type,
    bool preferSparseReducers
  ) {
    auto p = new ConcreteMonoLookup<BaseLookup, AllowRemovals, UseDivMask>
      (monoid, type, preferSparseReducers);
    return std::unique_ptr<MonoLookup>(p);
  }

  std::unique_ptr<MonoLookup> createGeneral(
    const MonoLookup::Monoid& monoid,
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

  class ConcreteFactory : public MonoLookup::Factory {
  public:
    ConcreteFactory(const Monoid& monoid, int type): 
      mMonoid(monoid),
      mType(type)
    {}

    virtual std::unique_ptr<MonoLookup> create(
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

std::unique_ptr<MonoLookup::Factory> MonoLookup::makeFactory(
  const PolyRing& ring,
  int type
) {
  return std::unique_ptr<Factory>(new ConcreteFactory(ring.monoid(), type));
}

void MonoLookup::displayMonoLookupTypes(std::ostream& out) {
  out << "Mono Lookup Types:" << std::endl;
  out << "  1   divlist+divmask" << std::endl;
  out << "  2   kdtree+divmask" << std::endl;
  out << "  3   divlist" << std::endl;
  out << "  4   kdtree" << std::endl;
}

MATHICGB_NAMESPACE_END
