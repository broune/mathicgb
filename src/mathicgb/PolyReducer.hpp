// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_POLY_REDUCER_GUARD
#define MATHICGB_POLY_REDUCER_GUARD

#include "TypicalReducer.hpp"

MATHICGB_NAMESPACE_BEGIN

class PolyReducer : public TypicalReducer {
public:
  PolyReducer(const PolyRing *R);

  virtual ~PolyReducer();

  virtual std::string description() const { return "poly reducer"; }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term& result);
  virtual void removeLeadTerm();

  void value(Poly &result); // keep extracting lead term until done
  void dump() const;

  virtual size_t getMemoryUse() const;

protected:
  void resetReducer();

private:
  void insert(const_term multiplier, Poly::const_iterator first, Poly::const_iterator last);

  const PolyRing *R;
  Poly *f;
  Poly::iterator f_iter;

  size_t mMemUsage;
};

MATHICGB_NAMESPACE_END
#endif
