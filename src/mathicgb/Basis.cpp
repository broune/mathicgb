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

std::unique_ptr<Basis> Basis::parse(std::istream& in)
{
  PolyRing *R = PolyRing::read(in); // todo: fix this leak
  size_t npolys;
  in >> npolys;
  auto result = make_unique<Basis>(*R);
  for (size_t j = 0; j < npolys; ++j) {
    auto g = make_unique<Poly>(*R);
    while (std::isspace(in.peek()))
      in.get();
    g->parse(in);
    result->insert(std::move(g));
  }
  return result;
}

void Basis::display(std::ostream& out, bool printComponent) const
{
  mRing.write(out);
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
