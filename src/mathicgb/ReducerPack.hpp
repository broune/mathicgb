// Copyright 2011 Bjarke Roune, Michael E. Stillman

#ifndef _reducer_pack_h_
#define _reducer_pack_h_

#include <memtailor.h>
#include <mathic.h>

#include "TypicalReducer.hpp"
#include "ReducerHelper.hpp"

template<template<typename> class Queue>
class ReducerPack : public TypicalReducer {
public:
  ReducerPack(const PolyRing& ring);
  virtual ~ReducerPack();

  virtual std::string description() const {
    return mQueue.getName() + "-packed";
  }

  virtual void insertTail(const_term multiplier, const Poly* f);
  virtual void insert(monomial multiplier, const Poly* f);

  virtual bool leadTerm(const_term& result);
  virtual void removeLeadTerm();

  virtual size_t getMemoryUse() const;

protected:
  virtual void resetReducer();

private:
  // Represents a term multiple of a polynomial, 
  // together with a current term of the multiple.
public:
  struct MultipleWithPos {
    MultipleWithPos(const Poly& poly, const_term multiple);

    Poly::const_iterator pos;
    Poly::const_iterator const end;
    const_term const multiple;

    // invariant: current is the monomial product of multiple.monom 
    // and pos.getMonomial().
    monomial current;

    // Ensures the invariant, so sets current to the product of
    // multiple.monom and pos.getMonomial().
    void computeCurrent(const PolyRing& ring);
    void currentCoefficient(const PolyRing& ring, coefficient& coeff);
    void destroy(const PolyRing& ring);
  };

  class Configuration : public ReducerHelper::PlainConfiguration {
  public:
    typedef MultipleWithPos* Entry;
    Configuration(const PolyRing& ring) : PlainConfiguration(ring) {}
    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monomialLT(a->current, b->current);
    }
  };
private:

  class MonomialFree;
  
  const PolyRing& mRing;
  term mLeadTerm;
  bool mLeadTermKnown;
  Queue<Configuration> mQueue;
  memt::BufferPool mPool;
};

extern int tracingLevel;

template<template<typename> class Q>
ReducerPack<Q>::ReducerPack(const PolyRing& ring):
  mRing(ring),
  mLeadTerm(0, mRing.allocMonomial()),
  mLeadTermKnown(false),
  mQueue(Configuration(ring)),
  mPool(sizeof(MultipleWithPos))
{
}

template<template<typename> class Q>
class ReducerPack<Q>::MonomialFree
{
public:
  MonomialFree(const PolyRing& ring): mRing(ring) {}

  bool proceed(MultipleWithPos* entry)
  {
    entry->destroy(mRing);
    return true;
  }
private:
  const PolyRing& mRing;
};

template<template<typename> class Q>
ReducerPack<Q>::~ReducerPack()
{
  resetReducer();
  mRing.freeMonomial(mLeadTerm.monom);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
template<template<typename> class Q>
void ReducerPack<Q>::insertTail(const_term multiple, const Poly* poly)
{
  if (poly->nTerms() <= 1)
    return;
  mLeadTermKnown = false;

  MultipleWithPos* entry =
    new (mPool.alloc()) MultipleWithPos(*poly, multiple);
  ++entry->pos;
  entry->computeCurrent(poly->ring());
  mQueue.push(entry);
}

template<template<typename> class Q>
void ReducerPack<Q>::insert(monomial multiple, const Poly* poly)
{
  if (poly->isZero())
    return;
  mLeadTermKnown = false;

  // todo: avoid multiplication by 1
  term termMultiple(1, multiple);
  MultipleWithPos* entry =
    new (mPool.alloc()) MultipleWithPos(*poly, termMultiple);
  entry->computeCurrent(poly->ring());
  mQueue.push(entry);
}

template<template<typename> class Q>
ReducerPack<Q>::MultipleWithPos::MultipleWithPos
(const Poly& poly, const_term multipleParam):
  pos(poly.begin()),
  end(poly.end()),
  multiple(ReducerHelper::allocTermCopy(poly.ring(), multipleParam)),
  current(poly.ring().allocMonomial()) {}

template<template<typename> class Q>
void ReducerPack<Q>::MultipleWithPos::computeCurrent(const PolyRing& ring) {
  ring.monomialMult(multiple.monom, pos.getMonomial(), current);  
}

template<template<typename> class Q>
void ReducerPack<Q>::MultipleWithPos::currentCoefficient
(const PolyRing& ring, coefficient& coeff) {
  ring.coefficientMult(multiple.coeff, pos.getCoefficient(), coeff);
}

template<template<typename> class Q>
void ReducerPack<Q>::MultipleWithPos::destroy(const PolyRing& ring) {
  ring.freeMonomial(current);
  ConstMonomial& monom = const_cast<ConstMonomial&>(multiple.monom);
  ring.freeMonomial(monom.castAwayConst());

  // Call the destructor to destruct the iterators into std::vector.
  // In debug mode MSVC puts those in a linked list and the destructor
  // has to be called since it takes an iterator off the list. We had
  // memory corruption problems before doing this.
  this->~MultipleWithPos();
}

template<template<typename> class Q>
bool ReducerPack<Q>::leadTerm(const_term& result)
{
  if (mLeadTermKnown) {
    result = mLeadTerm;
    return true;
  }

  do {
    if (mQueue.empty())
      return false;
    MultipleWithPos* entry = mQueue.top();
    mLeadTerm.monom.swap(entry->current);
    entry->currentCoefficient(mRing, mLeadTerm.coeff);
    
    while (true) {
      ++entry->pos;
      if (entry->pos == entry->end) {
        mQueue.pop();
        entry->destroy(mRing);
        mPool.free(entry);
      } else {
        entry->computeCurrent(mRing);
        mQueue.decreaseTop(entry);
      }
      
      if (mQueue.empty())
        break;
      
      entry = mQueue.top();
      if (!mRing.monomialEQ(entry->current, mLeadTerm.monom))
        break;
      coefficient coeff;
      entry->currentCoefficient(mRing, coeff);
      mRing.coefficientAddTo
        (mLeadTerm.coeff, const_cast<const coefficient&>(coeff));
    }
  } while (mRing.coefficientIsZero(mLeadTerm.coeff));

  result = mLeadTerm;
  mLeadTermKnown = true;
  return true;
}

template<template<typename> class Q>
void ReducerPack<Q>::removeLeadTerm()
{
  if (!mLeadTermKnown) {
    const_term dummy;
    leadTerm(dummy);
  }
  mLeadTermKnown = false;
}

template<template<typename> class Q>
void ReducerPack<Q>::resetReducer()
{
  MonomialFree freeer(mRing);
  mQueue.forAll(freeer);
  mQueue.clear();
}

template<template<typename> class Q>
size_t ReducerPack<Q>::getMemoryUse() const
{
  return
    TypicalReducer::getMemoryUse() +
    mQueue.getMemoryUse() +
    mPool.getMemoryUse();
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
