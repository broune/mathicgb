// Copyright 2011 Michael E. Stillman

#ifndef _polyHeap_h_
#define _polyHeap_h_

#include "TypicalReducer.hpp"

struct heap_term {
  monomial actual;  // multiplier * *first
  const_term multiplier;
  Poly::const_iterator first;
  Poly::const_iterator last;
};

class heapCompareFcn {
  const PolyRing *R;
public:
  heapCompareFcn(const PolyRing *R0) : R(R0) {}
  bool operator()(heap_term &a, heap_term &b);
};

class PolyHeap : public TypicalReducer {
public:
  PolyHeap(const PolyRing *R);
  ~PolyHeap() {}

  virtual std::string description() const { return "heap"; }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term& result);
  virtual void removeLeadTerm();

  void value(Poly &result); // keep extracting lead term until done
  void dump() const;

  static size_t stats_static_n_compares; // static so that heap compares can be easily counted

  void insert(const_term multiplier, Poly::const_iterator first, Poly::const_iterator last);

  size_t getMemoryUse() const;

protected:
  void resetReducer();

private:
  bool extractLeadTerm(const_term &result);
  // returns true if there is a term to extract

  const PolyRing *R;
  heapCompareFcn heapCompare;
  std::vector<heap_term> terms_;  // the actual heap

  bool lead_is_computed;
  const_term lead_term;  // if it has been computed, this will be set.
  // monomial points somewhere into 'monoms'.

  void extract1(const_term &t); // assumes: heap is not empty
  bool computeLeadTerm();
};

inline bool heapCompareFcn::operator()(heap_term &a, heap_term &b)
{
  PolyHeap::stats_static_n_compares++;
  return R->monomialLT(a.actual, b.actual);
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
