// Copyright 2011 Michael E. Stillman

#ifndef _polyReducer_h_
#define _polyReducer_h_

#include "Reducer.hpp"

class PolyReducer : public Reducer {
public:
  PolyReducer(const PolyRing *R);

  virtual ~PolyReducer();

  virtual std::string description() const { return "poly reducer"; }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  bool findLeadTerm(const_term &result);
  void removeLeadTerm();

  void value(Poly &result); // keep extracting lead term until done
  void dump() const;

  size_t getMemoryUse() const { return mMemUsage; } // this one reports high water usage
protected:
  void resetReducer();

private:
  void insert(const_term multiplier, Poly::const_iterator first, Poly::const_iterator last);

  const PolyRing *R;
  Poly *f;
  Poly::iterator f_iter;

  size_t mMemUsage;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
