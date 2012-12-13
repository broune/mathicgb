#ifndef MATHICGB_TYPICAL_REDUCER_GUARD
#define MATHICGB_TYPICAL_REDUCER_GUARD

#include "Reducer.hpp"
#include "Poly.hpp"
#include "PolyRing.hpp"
class GroebnerBasis;
class PolyBasis;

/** Uses the template method pattern (not C++ templates) to implement
  reduction with some steps left to sub-classes.

  The template method pattern defines an algorithm that calls several
  virtual methods. Sub-classes can then redefine those methods to
  change some parts of the algorithm without recoding the whole
  thing. The word "template" here has nothing to do with C++
  templates. See http://en.wikipedia.org/wiki/Template_method_pattern
*/
class TypicalReducer : public Reducer {
public:
  virtual Poly* regularReduce(
    const_monomial sig,
    const_monomial multiple,
    size_t basisElement,
    const GroebnerBasis& basis);

  virtual std::unique_ptr<Poly> classicReduce
  (const Poly& poly, const PolyBasis& basis);

  virtual std::unique_ptr<Poly> classicTailReduce
  (const Poly& poly, const PolyBasis& basis);

  virtual std::unique_ptr<Poly> classicReduceSPoly
    (const Poly& a, const Poly& b, const PolyBasis& basis);

  virtual void classicReduceSPolySet
  (std::vector<std::pair<size_t, size_t> >& spairs,
   const PolyBasis& basis,
   std::vector<std::unique_ptr<Poly> >& reducedOut);

  virtual void classicReducePolySet
  (const std::vector<std::unique_ptr<Poly> >& polys,
   const PolyBasis& basis,
   std::vector<std::unique_ptr<Poly> >& reducedOut);

  virtual void setMemoryQuantum(size_t quantum);

protected:
  // These are the methods that sub-classes define in order to carry
  // out sub-steps in the reduction.
  virtual void insertTail(const_term multiplier, const Poly *f) = 0;
  virtual void insert(monomial multiplier, const Poly *f) = 0;
  virtual bool leadTerm(const_term& result) = 0;
  virtual void removeLeadTerm() = 0;
  virtual void resetReducer() = 0;

  virtual size_t getMemoryUse() const;

  // Sub-classes can use this
  memt::Arena mArena;

private:
  void reset();
  std::unique_ptr<Poly> classicReduce(const PolyBasis& basis);
  std::unique_ptr<Poly> classicReduce
    (std::unique_ptr<Poly> partialResult, const PolyBasis& basis);
};

#endif
