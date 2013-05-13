// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_POLY_GEO_BUCKET_GUARD
#define MATHICGB_POLY_GEO_BUCKET_GUARD

#include "TypicalReducer.hpp"

MATHICGB_NAMESPACE_BEGIN

/// @todo: shouldn't be a macro
#define GEOHEAP_SIZE 15

class PolyGeoBucket : public TypicalReducer {
public:
  PolyGeoBucket(const PolyRing *R);
  ~PolyGeoBucket();

  virtual std::string description() const { return "geo buckets"; }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term& result);
  virtual void removeLeadTerm();

  void value(Poly &result); // keep extracting lead term until done
  void resetReducer();

  void dump() const; // Used for debugging

  size_t getMemoryUse() const;

private:
  void insert(const_term multiplier, Poly::iterator first, Poly::iterator last);

  struct heap_record {
    Poly *poly;
    Poly::iterator first;
    Poly::iterator last;
  };

  void add_to(heap_record &a, heap_record &b);
  // sets a to a+b, sets b to 0.

private:
  const PolyRing *R;
  heap_record heap[GEOHEAP_SIZE];
  int top_of_heap;
  int lead; // where the lead monomial is located (i.e. heap[lead]).  = -1 means not known.

  void insert(heap_record &a);
  bool computeLeadTerm();
  void update_size_stats();

};

MATHICGB_NAMESPACE_END
#endif
