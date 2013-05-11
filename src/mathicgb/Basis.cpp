#include "stdinc.h"
#include "Basis.hpp"

#include "PolyRing.hpp"
#include "Poly.hpp"
#include "MathicIO.hpp"
#include <ostream>
#include <istream>
#include <iostream>
#include <cctype>

void Basis::insert(std::unique_ptr<Poly>&& p) {
  MATHICGB_ASSERT(p.get() != 0);
  MATHICGB_ASSERT(p->termsAreInDescendingOrder());
  mGenerators.push_back(std::move(p));
}

void Basis::sort() {
  const auto& monoid = ring().monoid();
  const auto cmp = [&monoid](
    const std::unique_ptr<Poly>& a,
    const std::unique_ptr<Poly>& b
  ) {
    return monoid.lessThan(a->getLeadMonomial(), b->getLeadMonomial());
  };
  std::sort(mGenerators.begin(), mGenerators.end(), cmp);
}

void Basis::display(std::ostream& out, bool printComponent, bool componentIncreasingDesired) const
{
  mRing.write(out, componentIncreasingDesired);
  out << '\n' << mGenerators.size() << '\n';
  for (size_t i = 0; i < mGenerators.size(); ++i) {
    mGenerators[i]->display(out, printComponent);
    out << '\n';
  }
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
