// Copyright 2011 Michael E. Stillman

#ifndef _monomial_module_h_
#define _monomial_module_h_

#include "PolyRing.hpp"
#include <vector>
#include <iostream>
#include <algorithm>

class MonomialTableArray
{
public:
  typedef size_t V; // value type

  static MonomialTableArray *makeMonomialTableArray(int type, // 0=naive, 1=DivList, 2=KDTree
                                                    const PolyRing *R,
                                                    // for each type, only some of the following are relevant
                                                    bool use_cache,
                                                    bool remin);

  MonomialTableArray(const PolyRing *R0) : R(R0) {}
  virtual ~MonomialTableArray() {};

  virtual bool insert(const_monomial m, V val) = 0;
  // returns true if the monomial actually needs to be inserted.
  // If the monomial is inserted, the caller agrees to keep that monomial
  // active until it is removed from the table.
  //TODO: At some point: deal with removals too

  // Adds a new component, increasing the number of monomial tables
  // in the array by one. Its index is one more than the previous
  // maximum index.
  virtual void addComponent() = 0;

  bool member(const_monomial m) {
    size_t dummy;
    return member(m, dummy);
  }
  virtual bool member(const_monomial m, V &result) = 0;

  struct Stats {
    size_t n_actual_inserts;
    size_t n_calls_member;
    size_t n_calls_insert;
    size_t n_compares;

    // ratio of non-zero exponents. 0.75 means a
    // quarter of all exponents are zero.
    double denseness;

    Stats() : n_actual_inserts(0),
              n_calls_member(0),
              n_calls_insert(0),
              n_compares(0),
              denseness(0.0) {}
  };

  virtual void getStats(Stats &stats) const = 0;

  virtual void display(std::ostream &o, int level) const = 0;

  virtual void printFrobbyM2Format
    (std::ostream& out, size_t component) const = 0;

  virtual void displayStats(std::ostream &o) const = 0;

  virtual size_t n_elems() const = 0;

  void dump(int level=0) const
  {
    // display on stderr the table.
    displayStats(std::cerr);
    if (level > 0) display(std::cerr, level-1);
  }

  virtual std::string description() const = 0;

  virtual size_t getMemoryUse() const = 0;

  // Choosing one
  static int displayMTTypes(std::ostream &o); // returns n s.t. 0..n-1 are valid types
  static std::unique_ptr<MonomialTableArray>
    make(const PolyRing *R, int typ, size_t components, bool allowRemovals);

protected:
  class MonomialCompare {
  public:
    MonomialCompare(const PolyRing& ring): mRing(ring) {}
    bool operator()(const_monomial a, const_monomial b) {
      return mRing.monomialLT(a, b);
    }
  private:
    const PolyRing& mRing;
  };

  const PolyRing *R;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
