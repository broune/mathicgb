// Copyright 2011 Michael E. Stillman
#include "stdinc.h"
#include "PolyRing.hpp"

#include "Poly.hpp"
#include "FreeModuleOrder.hpp"
#include "Ideal.hpp"
#include <ostream>
#include <istream>
#include <iostream>
#include <cctype>

Ideal::~Ideal()
{
  for (size_t i = 0; i<mGenerators.size(); i++)
    delete mGenerators[i];
}

void Ideal::insert(std::unique_ptr<Poly> p) {
  MATHICGB_ASSERT(p.get() != 0);
  MATHICGB_ASSERT(p->termsAreInDescendingOrder());
  mGenerators.reserve(mGenerators.size() + 1);
  mGenerators.push_back(p.release());
}

namespace {
  class IdealSort {
  public:
    IdealSort(const FreeModuleOrder& order): mOrder(order) {}
    bool operator()(const Poly* a, const Poly* b) {
      return mOrder.signatureCompare
        (a->getLeadMonomial(), b->getLeadMonomial()) == LT;
    }

  private:
    const FreeModuleOrder& mOrder;
  };
}

void Ideal::sort(FreeModuleOrder& order) {
  IdealSort cmp(order);
  std::sort(mGenerators.begin(), mGenerators.end(), cmp);
}

std::unique_ptr<Ideal> Ideal::parse(std::istream& in)
{
  PolyRing *R = PolyRing::read(in); // todo: fix this leak
  size_t npolys;
  in >> npolys;
  auto result = make_unique<Ideal>(*R);
  for (size_t j = 0; j < npolys; ++j) {
    auto g = make_unique<Poly>(*R);
    while (std::isspace(in.peek()))
      in.get();
    g->parse(in);
    result->insert(std::move(g));
  }
  return result;
}

void Ideal::display(std::ostream &o, bool print_comp) const
{
  mRing.write(o);
  o << std::endl;
  o << mGenerators.size() << std::endl;
  for (size_t i = 0; i<mGenerators.size(); i++)
    {
      mGenerators[i]->display(o, print_comp);
      o << std::endl;
    }
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
