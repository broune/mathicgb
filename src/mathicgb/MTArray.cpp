// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "MTArray.hpp"

#include "StaticMonoLookup.hpp"
#include "MathicIO.hpp"

MATHICGB_NAMESPACE_BEGIN

template<
  bool UseKDTree,
  bool AllowRemovals,
  bool UseDivMask
>
class ConcreteModuleMonoSet : public MonomialTableArray {
public:
  typedef PolyRing::Monoid Monoid;
  struct NoData {};
  typedef StaticMonoLookup<UseKDTree, NoData, AllowRemovals, UseDivMask>
    Lookup;
  typedef typename Lookup::Entry Entry;

  ConcreteModuleMonoSet(const Monoid& monoid, const size_t componentCount):
    mMonoid(monoid),
    mComponentCount(componentCount),
    mLookups(reinterpret_cast<Lookup*>(
      new char[sizeof(Lookup) * componentCount]
    ))
  {
    for (size_t component = 0; component < componentCount; ++component)
      new (&mLookups[component]) Lookup(monoid);
  }

  virtual ~ConcreteModuleMonoSet() {
    delete[] reinterpret_cast<char*>(mLookups);
  }

  virtual bool insert(const_monomial m) {
    const auto c = monoid().component(m);
    MATHICGB_ASSERT(c < componentCount());
    if (member(m))
      return false;
    if (AllowRemovals)
      mLookups[c].removeMultiples(m);
    mLookups[c].insert(m, NoData());
    return true;
  }

  virtual bool member(const_monomial m) {
    const auto c = monoid().component(m);
    MATHICGB_ASSERT(c < componentCount());
    return mLookups[c].divisor(m) != 0;
  }

  virtual std::string name() const {
    return Lookup(monoid()).getName();
  }

  virtual void display(std::ostream& out) const {
    std::vector<const_monomial> monomials;
    for (size_t c = 0; c < componentCount(); ++c) {
      const Lookup& lookup = mLookups[c];
      if (lookup.size() == 0)
        continue;
      out << "  " << c << ": ";
      monomials.clear();
      lookup.forAll([&](const Entry& entry) {
        monomials.emplace_back(Monoid::toOld(entry.mono()));
      });
      const auto& monoid = this->monoid(); // workaround for gcc 4.5.3 issue
      const auto cmp = [&](const_monomial a, const_monomial b) {
        return monoid.lessThan(a, b);
      };
      std::sort(monomials.begin(), monomials.end(), cmp);
      for (auto mono = monomials.cbegin(); mono != monomials.cend(); ++mono) {
        MathicIO<>().writeMonomial(monoid, false, *mono, out);
        out << "  ";
      }
      out << '\n';
    }
  }

  virtual void getMonomials(std::vector<const_monomial>& monomials) const {
    for (size_t c = 0; c < componentCount(); ++c) {
      mLookups[c].forAll(
        [&](const Entry& entry) {
          monomials.push_back(Monoid::toOld(entry.mono()));
        }
      );
    }
  }
  
  virtual size_t n_elems() const {
    size_t count = 0;
    for (size_t c = 0; c < componentCount(); ++c)
      count += mLookups[c].size();
    return count;
  }

  virtual size_t getMemoryUse() const {
    size_t count = 0;
    for (size_t c = 0; c < componentCount(); ++c)
      count += mLookups[c].getMemoryUse();
    return count;
  }

  const Monoid& monoid() const {return mMonoid;}
  size_t componentCount() const {return mComponentCount;}

private:
  const Monoid& mMonoid;

  // We cannot use std::vector as Lookup might not be copyable or movable.
  const size_t mComponentCount;
  Lookup* const mLookups;
};

int MonomialTableArray::displayMTTypes(std::ostream &o)
 // returns n s.t. 0..n-1 are valid types
{
  o << "Monomial table types:" << std::endl;
  o << "  1   divlist+divmask" << std::endl;
  o << "  2   kdtree+divmask" << std::endl;
  o << "  3   divlist" << std::endl;
  o << "  4   kdtree" << std::endl;
  return 5;
}

namespace {
  class ModuleMonoSetFactory {
  public:
    typedef PolyRing::Monoid Monoid;

    std::unique_ptr<MonomialTableArray> create(
      const Monoid& monoid,
      int type,
      size_t componentCount,
      bool allowRemovals
    ) {
      const Params params = {monoid, type, componentCount};
      return staticMonoLookupCreate
        <Create, std::unique_ptr<MonomialTableArray>>
        (type, allowRemovals, params);
    }

  private:
    struct Params {
      const Monoid& monoid;
      int type;
      size_t componentCount;
    };

    template<bool UseKDTree, bool AllowRemovals, bool UseDivMask>
    struct Create {
      static std::unique_ptr<MonomialTableArray> create(const Params& params) {
        return make_unique
          <ConcreteModuleMonoSet<UseKDTree, AllowRemovals, UseDivMask>>
          (params.monoid, params.componentCount);
      }
    };
  };
}

std::unique_ptr<MonomialTableArray> MonomialTableArray::make
  (const PolyRing *R, int typ, size_t components, bool allowRemovals)
{
  return ModuleMonoSetFactory().create(R->monoid(), typ, components, allowRemovals);
}

MATHICGB_NAMESPACE_END
