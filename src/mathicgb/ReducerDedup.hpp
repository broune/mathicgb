// Copyright 2011 Bjarke Roune, Michael E. Stillman

#ifndef _reducer_dedup_h_
#define _reducer_dedup_h_

#include <memtailor.h>
#include <mathic.h>

#include "Reducer.hpp"
#include "ReducerHelper.hpp"

template<template<typename ConfigType> class Queue> class ReducerDedup;

template<template<typename> class Queue>
class ReducerDedup : public Reducer {
public:
  ReducerDedup(const PolyRing& R);
  virtual ~ReducerDedup();

  virtual std::string description() const { 
    return mQueue.getName() + "-dedup"; 
  }

  virtual void insertTail(const_term multiplier, const Poly *f);
  virtual void insert(monomial multiplier, const Poly *f);

  virtual bool findLeadTerm(const_term &result);
  virtual void removeLeadTerm();

  virtual size_t getMemoryUse() const;

protected:
  virtual void resetReducer();

public:
  // This Configuration is designed to work with
  // mathic::TourTree, mathic::Heap, and mathic::Geobucket

  class Configuration : public ReducerHelper::DedupConfiguration {
  public:
    typedef term Entry;
    Configuration(const PolyRing& ring): DedupConfiguration(ring) {}
    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monomialCompare(a.monom, b.monom);
    }
    Entry deduplicate(Entry a, Entry b) const {
      // change a.coeff, and free b.monom
      ring().coefficientAddTo(a.coeff, b.coeff);
      ring().freeMonomial(b.monom);
      return a;
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
ReducerDedup<Q>::ReducerDedup(const PolyRing& ring):
  Reducer(),
  mRing(ring),
  mLeadTerm(0, mRing.allocMonomial()),
  mLeadTermKnown(false),
  mQueue(Configuration(ring))
{
}

template<template<typename> class Q>
class ReducerDedup<Q>::MonomialFree
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
ReducerDedup<Q>::~ReducerDedup()
{
  resetReducer();
  mRing.freeMonomial(mLeadTerm.monom);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
template<template<typename> class Q>
void ReducerDedup<Q>::insertTail(const_term multiple, const Poly* poly)
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
void ReducerDedup<Q>::insert(monomial multiple, const Poly* poly)
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
bool ReducerDedup<Q>::findLeadTerm(const_term& result)
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
void ReducerDedup<Q>::removeLeadTerm()
{
  if (!mLeadTermKnown) {
    const_term dummy;
    findLeadTerm(dummy);
  }
  mLeadTermKnown = false;
}

template<template<typename> class Q>
void ReducerDedup<Q>::resetReducer()
{
  MonomialFree freeer(mRing);
  //mQueue.forAll(freeer);
  //  mQueue.clear();
}

template<template<typename> class Q>
size_t ReducerDedup<Q>::getMemoryUse() const
{
  return mQueue.getMemoryUse();
}


#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
