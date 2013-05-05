#include "stdinc.h"
#include "Basis.hpp"

#include "PolyRing.hpp"
#include "Poly.hpp"
#include "FreeModuleOrder.hpp"
#include <ostream>
#include <istream>
#include <iostream>
#include <cctype>

Basis::~Basis()
{
  for (size_t i = 0; i<mGenerators.size(); i++)
    delete mGenerators[i];
}

void Basis::insert(std::unique_ptr<Poly> p) {
  MATHICGB_ASSERT(p.get() != 0);
  MATHICGB_ASSERT(p->termsAreInDescendingOrder());
  mGenerators.reserve(mGenerators.size() + 1);
  mGenerators.push_back(p.release());
}

namespace {
  class BasisSort {
  public:
    BasisSort(const FreeModuleOrder& order): mOrder(order) {}
    bool operator()(const Poly* a, const Poly* b) {
      return mOrder.signatureCompare
        (a->getLeadMonomial(), b->getLeadMonomial()) == LT;
    }

  private:
    const FreeModuleOrder& mOrder;
  };
}

void Basis::sort(FreeModuleOrder& order) {
  BasisSort cmp(order);
  std::sort(mGenerators.begin(), mGenerators.end(), cmp);
}

auto Basis::parse(std::istream& in) -> Parsed
{
  auto r = PolyRing::read(in);
  auto ring = make_unique<PolyRing>(std::move(*r.first));
  delete r.first;

  auto basis = make_unique<Basis>(*ring);
  auto processor =
    make_unique<MonoProcessor<Monoid>>(ring->monoid(), r.second.first, false);

  size_t polyCount;
  in >> polyCount;
  for (size_t i = 0; i < polyCount; ++i) {
    auto poly = make_unique<Poly>(*ring);
    while (std::isspace(in.peek()))
      in.get();
    poly->parse(in);
    basis->insert(std::move(poly));
  }

  if (r.second.second) {
    Monoid::MonoVector schreyer(ring->monoid());
    for (size_t gen = 0; gen < basis->size(); ++gen)
      schreyer.push_back(basis->getPoly(gen)->getLeadMonomial());
    processor->setModuleAdjustments(std::move(schreyer));
  }

  return std::make_tuple(
    std::move(ring),
    std::move(basis),
    std::move(processor)
  );
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
