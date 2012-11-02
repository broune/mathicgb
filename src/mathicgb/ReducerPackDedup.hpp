// Copyright 2011 Bjarke Roune, Michael E. Stillman

#ifndef _reducer_pack_dedup_h_
#define _reducer_pack_dedup_h_

#include <memtailor.h>
#include <mathic.h>

#include "TypicalReducer.hpp"
#include "ReducerHelper.hpp"

template<template<typename> class Queue>
class ReducerPackDedup : public TypicalReducer {
public:
  ReducerPackDedup(const PolyRing& ring);
  virtual ~ReducerPackDedup();

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
    void addCurrentCoefficient(const PolyRing& ring, coefficient& coeff);
    void destroy(const PolyRing& ring);

    // Points to a circular list of entries that have the same current
    // monomial. If no other such entries have been discovered, then
    // chain points to this object itself. We use a circular linked list
    // as those allow merging in O(1) time.
    MultipleWithPos* chain;
    void mergeChains(MultipleWithPos& entry) {
      // This only works if *this and entry are not already in the
      // same chain!
      std::swap(chain, entry.chain);
    }
  };

  class Configuration : public ReducerHelper::DedupConfiguration {
  public:
    typedef MultipleWithPos* Entry;
    Configuration(const PolyRing& ring): DedupConfiguration(ring) {}
    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monomialCompare(a->current, b->current);
    }
    Entry deduplicate(Entry a, Entry b) const {
      a->mergeChains(*b);
      return a;
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

template<template<typename> class Q>
ReducerPackDedup<Q>::ReducerPackDedup(const PolyRing& ring):
  mRing(ring),
  mLeadTerm(0, mRing.allocMonomial()),
  mLeadTermKnown(false),
  mQueue(Configuration(ring)),
  mPool(sizeof(MultipleWithPos)) {
}

template<template<typename> class Q>
class ReducerPackDedup<Q>::MonomialFree
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
ReducerPackDedup<Q>::~ReducerPackDedup()
{
  resetReducer();
  mRing.freeMonomial(mLeadTerm.monom);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
template<template<typename> class Q>
void ReducerPackDedup<Q>::insertTail(const_term multiple, const Poly* poly)
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
void ReducerPackDedup<Q>::insert(monomial multiple, const Poly* poly)
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
ReducerPackDedup<Q>::MultipleWithPos::MultipleWithPos
(const Poly& poly, const_term multipleParam):
  pos(poly.begin()),
  end(poly.end()),
  multiple(ReducerHelper::allocTermCopy(poly.ring(), multipleParam)),
  current(poly.ring().allocMonomial()),
  chain(this) {}

template<template<typename> class Q>
void ReducerPackDedup<Q>::MultipleWithPos::computeCurrent(const PolyRing& ring) {
  ring.monomialMult(multiple.monom, pos.getMonomial(), current);  
}

template<template<typename> class Q>
void ReducerPackDedup<Q>::MultipleWithPos::currentCoefficient
(const PolyRing& ring, coefficient& coeff) {
  ring.coefficientMult(multiple.coeff, pos.getCoefficient(), coeff);
}

template<template<typename> class Q>
void ReducerPackDedup<Q>::MultipleWithPos::addCurrentCoefficient
(const PolyRing& ring, coefficient& coeff) {
  coefficient tmp;
  ring.coefficientMult(multiple.coeff, pos.getCoefficient(), tmp);
  ring.coefficientAddTo(coeff, tmp);
}

template<template<typename> class Q>
void ReducerPackDedup<Q>::MultipleWithPos::destroy(const PolyRing& ring) {
  MultipleWithPos* entry = this;
  do {
    ring.freeMonomial(entry->current);
    ConstMonomial& monom = const_cast<ConstMonomial&>(entry->multiple.monom);
    ring.freeMonomial(monom.castAwayConst());
    MultipleWithPos* next = entry->chain;
    MATHICGB_ASSERT(next != 0);

    // Call the destructor to destruct the iterators into std::vector.
    // In debug mode MSVC puts those in a linked list and the destructor
    // has to be called since it takes an iterator off the list. We had
    // memory corruption problems before doing this.
    entry->~MultipleWithPos();

    entry = next;
  } while (entry != this);
}

template<template<typename> class Q>
bool ReducerPackDedup<Q>::leadTerm(const_term& result)
{
  if (mLeadTermKnown) {
    result = mLeadTerm;
    return true;
  }

  do {
    if (mQueue.empty())
      return false;
    MultipleWithPos* entry = mQueue.top();
    entry->currentCoefficient(mRing, mLeadTerm.coeff);
    while (true) {
      // store the chained elements
      MultipleWithPos* const chainBegin = entry->chain;
      MultipleWithPos* const chainEnd = entry; // the list is circular
      entry->chain = entry; // detach any chained elements

      // handle the entry itself
      mLeadTerm.monom.swap(entry->current);
      ++entry->pos;
      if (entry->pos == entry->end) {
        mQueue.pop();
        entry->destroy(mRing);
        mPool.free(entry);
      } else {
        entry->computeCurrent(mRing);
        // Inserted spans must be in descending order
        MATHICGB_ASSERT(mQueue.getConfiguration().ring().
          monomialLT(entry->current, mLeadTerm.monom));
        mQueue.decreaseTop(entry);
      }

      // handle any chained elements
      MultipleWithPos* chain = chainBegin;
      while (chain != chainEnd) {
        MATHICGB_ASSERT(chain != 0);
        MATHICGB_ASSERT(mRing.monomialEQ(chain->current, mLeadTerm.monom));

        MultipleWithPos* const next = chain->chain;
        chain->chain = chain; // detach from remaining chained elements

        chain->addCurrentCoefficient(mRing, mLeadTerm.coeff);
        ++chain->pos;
        if (chain->pos == chain->end) {
          chain->destroy(mRing);
          mPool.free(chain);
        } else {
          chain->computeCurrent(mRing);
          // Inserted spans must be in descending order
          MATHICGB_ASSERT(mQueue.getConfiguration().ring().
            monomialLT(chain->current, mLeadTerm.monom));
          mQueue.push(chain);
        }
        chain = next;
      }
      
      if (mQueue.empty())
        break;
      
      entry = mQueue.top();
      if (!mRing.monomialEQ(entry->current, mLeadTerm.monom))
        break;
      entry->addCurrentCoefficient(mRing, mLeadTerm.coeff);
    }
  } while (mRing.coefficientIsZero(mLeadTerm.coeff));

  result = mLeadTerm;
  mLeadTermKnown = true;
  return true;
}

template<template<typename> class Q>
void ReducerPackDedup<Q>::removeLeadTerm()
{
  if (!mLeadTermKnown) {
    const_term dummy;
    leadTerm(dummy);
  }
  mLeadTermKnown = false;
}

template<template<typename> class Q>
void ReducerPackDedup<Q>::resetReducer()
{
  MonomialFree freeer(mRing);
  mQueue.forAll(freeer);
  mQueue.clear();
}

template<template<typename> class Q>
size_t ReducerPackDedup<Q>::getMemoryUse() const
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

