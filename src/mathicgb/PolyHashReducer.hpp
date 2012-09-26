// Copyright 2011 Michael E. Stillman

#ifndef _polyHashReducer_h_
#define _polyHashReducer_h_

#include "TypicalReducer.hpp"
#include "PolyHashTable.hpp"

class PolyHashReducer : public TypicalReducer {
public:
  PolyHashReducer(const PolyRing *R);

  virtual ~PolyHashReducer();

  virtual std::string description() const { return "poly hash reducer"; }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term &result);
  virtual void removeLeadTerm();

  void value(Poly &result); // keep extracting lead term until done
  void dump() const;

  virtual size_t getMemoryUse() const;

protected:
  virtual void resetReducer();

private:
  void insert(const_term multiplier, Poly::const_iterator first, Poly::const_iterator last);

  typedef PolyHashTable::MonomialArray HashPoly;

  void merge(
    const HashPoly::const_iterator &fbegin,
    const HashPoly::const_iterator &fend,
    const HashPoly::const_iterator &gbegin,
    const HashPoly::const_iterator &gend,
    HashPoly &result,
    size_t &n_compares);

  const PolyRing* R_;
  HashPoly* f_;
  HashPoly::iterator f_iter_;

  PolyHashTable H_;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
