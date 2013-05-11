#ifndef MATHICGB_BASIS_GUARD
#define MATHICGB_BASIS_GUARD

#include "Poly.hpp"
#include "PolyRing.hpp"
#include <tuple>
#include <memory>
#include <algorithm>
#include <vector>

class Poly;

// Really: a list of polynomials
// BUT ALSO maybe: includes memory areas for the polynomials?
class Basis {
public:
  Basis(const PolyRing &R) : mRing(R) {}
  Basis(Basis&& basis): mRing(basis.ring()), mGenerators(std::move(basis.mGenerators)) {}

  void insert(std::unique_ptr<Poly>&& p);

  /// inverse operation to parse().
  void display(std::ostream &o, bool print_comp, bool componentIncreasingDesired) const;

  const PolyRing& ring() const { return mRing; }

  const PolyRing *getPolyRing() const { return &mRing; }
  const std::vector<std::unique_ptr<Poly>>& viewGenerators() {
    return mGenerators;
  }
  const Poly *getPoly(size_t i) const {
    MATHICGB_ASSERT(i < size());
    return mGenerators[i].get();
  }
  size_t size() const {return mGenerators.size();}
  bool empty() const {return mGenerators.empty();}
  void reserve(size_t size) {mGenerators.reserve(size);}

  void sort();

private:
  Basis(const Basis&); // not available

  const PolyRing& mRing;
  std::vector<std::unique_ptr<Poly>> mGenerators;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
