// Copyright 2011 Michael E. Stillman

#ifndef _ideal_h_
#define _ideal_h_

#include <memory>
#include <algorithm>
#include "PolyRing.hpp"

class Poly;
class FreeModuleOrder;

class Ideal {
  // Really: a list of polynomials
  // BUT ALSO maybe: includes memory areas for the polynomials?
public:
  Ideal(const PolyRing &R) : mRing(R) {}
  ~Ideal();

  void insert(std::unique_ptr<Poly> p);

  /// reads ring, #gens, each generator in turn
  static std::unique_ptr<Ideal> parse(std::istream &i);

  /// inverse operation to parse().
  void display(std::ostream &o, bool print_comp) const;

  const PolyRing& ring() const { return mRing; }

  const PolyRing *getPolyRing() const { return &mRing; }
  const std::vector<Poly *>& viewGenerators() { return mGenerators; }
  const Poly *getPoly(size_t i) const { return mGenerators[i]; }
  size_t size() const { return mGenerators.size(); }

  void sort(FreeModuleOrder& order);

private:
  const PolyRing& mRing;
  std::vector<Poly*> mGenerators;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
