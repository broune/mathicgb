// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "ReducerNoDedup.hpp"

#include "TypicalReducer.hpp"
#include "ReducerHelper.hpp"
#include <memtailor.h>
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

void reducerNoDedupDependency() {}

template<template<typename ConfigType> class Queue> class ReducerNoDedup;

template<template<typename> class Queue>
class ReducerNoDedup : public TypicalReducer {
public:
  ReducerNoDedup(const PolyRing& R);
  virtual ~ReducerNoDedup();

  virtual std::string description() const { 
    return mQueue.getName() + "-nodedup"; 
  }

  virtual void insertTail(const_term multiplier, const Poly *f);
  virtual void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term &result);
  virtual void removeLeadTerm();

  virtual size_t getMemoryUse() const;

protected:
  virtual void resetReducer();

public:
  // This Configuration is designed to work with
  // mathic::TourTree, mathic::Heap, and mathic::Geobucket

  class Configuration : public ReducerHelper::PlainConfiguration {
  public:
    typedef term Entry;
    Configuration(const PolyRing& ring): PlainConfiguration(ring) {}
    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monomialLT(a.monom, b.monom);
    }
  };

private:
  class MonomialFree;
  
  const PolyRing& mRing;
  term mLeadTerm;
  bool mLeadTermKnown;
  Queue<Configuration> mQueue;
};

template<template<typename> class Q>
ReducerNoDedup<Q>::ReducerNoDedup(const PolyRing& ring):
  mRing(ring),
  mLeadTerm(0, mRing.allocMonomial()),
  mLeadTermKnown(false),
  mQueue(Configuration(ring)) {
}

template<template<typename> class Q>
class ReducerNoDedup<Q>::MonomialFree
{
public:
  MonomialFree(const PolyRing& ring): mRing(ring) {}

  bool proceed(term entry)
  {
    mRing.freeMonomial(entry.monom);
    return true;
  }
private:
  const PolyRing& mRing;
};

template<template<typename> class Q>
ReducerNoDedup<Q>::~ReducerNoDedup()
{
  resetReducer();
  mRing.freeMonomial(mLeadTerm.monom);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
template<template<typename> class Q>
void ReducerNoDedup<Q>::insertTail(const_term multiple, const Poly* poly)
{
  if (poly->nTerms() <= 1)
    return;
  mLeadTermKnown = false;

  Poly::const_iterator i = poly->begin();
  for (++i; i != poly->end(); ++i)
    {
      term t;
      t.monom = mRing.allocMonomial();
      mRing.monomialMult(multiple.monom, i.getMonomial(), t.monom);
      mRing.coefficientMult(multiple.coeff, i.getCoefficient(), t.coeff);
      mQueue.push(t);
    }
}

template<template<typename> class Q>
void ReducerNoDedup<Q>::insert(monomial multiple, const Poly* poly)
{
  if (poly->isZero())
    return;
  mLeadTermKnown = false;

  for (Poly::const_iterator i = poly->begin(); i != poly->end(); ++i) {
    term t(i.getCoefficient(), mRing.allocMonomial());
    mRing.monomialMult(multiple, i.getMonomial(), t.monom);
    mQueue.push(t);
  }
}

template<template<typename> class Q>
bool ReducerNoDedup<Q>::leadTerm(const_term& result)
{
  if (mLeadTermKnown) {
    result = mLeadTerm;
    return true;
  }

  do {
    if (mQueue.empty())
      return false;
    mLeadTerm = mQueue.top();
    mQueue.pop();
    
    while (true) {
      if (mQueue.empty())
        break;
      
      term entry = mQueue.top();
      if (!mRing.monomialEQ(entry.monom, mLeadTerm.monom))
        break;
      mRing.coefficientAddTo(mLeadTerm.coeff, entry.coeff);
      mRing.freeMonomial(entry.monom);
      mQueue.pop();
    }
  } while (mRing.coefficientIsZero(mLeadTerm.coeff));

  result = mLeadTerm;
  mLeadTermKnown = true;
  return true;
}

template<template<typename> class Q>
void ReducerNoDedup<Q>::removeLeadTerm()
{
  if (!mLeadTermKnown) {
    const_term dummy;
    leadTerm(dummy);
  }
  mLeadTermKnown = false;
}

template<template<typename> class Q>
void ReducerNoDedup<Q>::resetReducer()
{
  MonomialFree freeer(mRing);
  //mQueue.forAll(freeer);
  //  mQueue.clear();
}

template<template<typename> class Q>
size_t ReducerNoDedup<Q>::getMemoryUse() const
{
  return TypicalReducer::getMemoryUse() + mQueue.getMemoryUse();
}

Reducer::Registration r1(
  "TourNoDedup",
  Reducer::Reducer_TourTree_NoDedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerNoDedup<mic::TourTree>>(ring);
  }
);

Reducer::Registration r7(
  "HeapNoDedup",
  Reducer::Reducer_Heap_NoDedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerNoDedup<mic::Heap>>(ring);
  }
);

Reducer::Registration r13(
  "GeoNoDedup",
  Reducer::Reducer_Geobucket_NoDedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerNoDedup<mic::Geobucket>>(ring);
  }
);


MATHICGB_NAMESPACE_END
