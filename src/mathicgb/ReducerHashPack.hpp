// Copyright 2011 Bjarke Roune, Michael E. Stillman

#ifndef _reducer_hash_pack_h
#define _reducer_hash_pack_h

#include <mathic.h>
#include <memtailor.h>

#include "TypicalReducer.hpp"
#include "ReducerHelper.hpp"
#include "PolyHashTable.hpp"

template<template<typename ConfigType> class Queue> class ReducerHashPack;

template<template<typename> class Queue>
class ReducerHashPack : public TypicalReducer {
public:
  ReducerHashPack(const PolyRing& R);
  virtual ~ReducerHashPack();

  virtual std::string description() const { 
    return mQueue.getName() + "-hashed-packed";
  }

  virtual void insertTail(const_term multiplier, const Poly *f);
  virtual void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term &result);
  virtual void removeLeadTerm();

  virtual size_t getMemoryUse() const;

protected:
  virtual void resetReducer();

private:
  // Represents a term multiple of a polynomial, 
  // together with a current term of the multiple.
  struct MultipleWithPos {
    MultipleWithPos(const Poly& poly, const_term multiple);

    Poly::const_iterator pos;
    Poly::const_iterator const end;
    const_term const multiple;
    monomial current; // multiple.monom * pos.getMonomial()
    PolyHashTable::node* node;

    void computeCurrent(const PolyRing& ring, monomial current);
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

  class MonomialFree;

  void insertEntry(MultipleWithPos* entry);

  const PolyRing& mRing;
  term mLeadTerm;
  bool mLeadTermKnown;
  Queue<Configuration> mQueue;
  PolyHashTable mHashTable;
  memt::BufferPool mPool;
};

template<template<typename> class Q>
ReducerHashPack<Q>::ReducerHashPack(const PolyRing& ring):
  mRing(ring),
  mLeadTerm(0, mRing.allocMonomial()),
  mLeadTermKnown(false),
  mQueue(Configuration(ring)),
  mHashTable(&ring, 10),
  mPool(sizeof(MultipleWithPos))
{
}

template<template<typename> class Q>
class ReducerHashPack<Q>::MonomialFree
{
public:
  MonomialFree(const PolyRing& ring): mRing(ring) {}

  bool proceed(MultipleWithPos* entry) {
    entry->destroy(mRing);
    return true;
  }
private:
  const PolyRing& mRing;
};

template<template<typename> class Q>
ReducerHashPack<Q>::~ReducerHashPack()
{
  resetReducer();
  mRing.freeMonomial(mLeadTerm.monom);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
template<template<typename> class Q>
void ReducerHashPack<Q>::insertTail(const_term multiple, const Poly* poly)
{
  MATHICGB_ASSERT(poly != 0);
  MATHICGB_ASSERT(&poly->ring() == &mRing);
  if (poly->nTerms() < 2)
    return;
  MultipleWithPos* entry =
    new (mPool.alloc()) MultipleWithPos(*poly, multiple);
  ++entry->pos;
  insertEntry(entry);
}

template<template<typename> class Q>
void ReducerHashPack<Q>::insert(monomial multiple, const Poly* poly)
{
  MATHICGB_ASSERT(poly != 0);
  MATHICGB_ASSERT(&poly->ring() == &mRing);
  if (poly->isZero())
    return;
  term termMultiple(1, multiple);
  insertEntry(new (mPool.alloc()) MultipleWithPos(*poly, termMultiple));
}

namespace {
  const_term allocTerm(const PolyRing& ring, const_term term) {
    monomial mono = ring.allocMonomial();
    ring.monomialCopy(term.monom, mono);
    return const_term(term.coeff, mono);
  }
}

template<template<typename> class Q>
ReducerHashPack<Q>::MultipleWithPos::MultipleWithPos
(const Poly& poly, const_term multiple):
  pos(poly.begin()),
  end(poly.end()),
  multiple(allocTerm(poly.ring(), multiple)),
  current(poly.ring().allocMonomial()),
  node(0) {}

template<template<typename> class Q>
void ReducerHashPack<Q>::MultipleWithPos::
computeCurrent(const PolyRing& ring, monomial current) {
  ring.monomialMult(multiple.monom, pos.getMonomial(), current);  
}

template<template<typename> class Q>
void ReducerHashPack<Q>::MultipleWithPos::currentCoefficient
(const PolyRing& ring, coefficient& coeff) {
  ring.coefficientMult(multiple.coeff, pos.getCoefficient(), coeff);
}

template<template<typename> class Q>
void ReducerHashPack<Q>::MultipleWithPos::destroy(const PolyRing& ring) {
  ring.freeMonomial(current);
  ring.freeMonomial(const_cast<ConstMonomial&>(multiple.monom).castAwayConst());

  // Call the destructor to destruct the iterators into std::vector.
  // In debug mode MSVC puts those in a linked list and the destructor
  // has to be called since it takes an iterator off the list. We had
  // memory corruption problems before doing this.
  this->~MultipleWithPos();
}

template<template<typename> class Q>
bool ReducerHashPack<Q>::leadTerm(const_term& result)
{
  if (mLeadTermKnown) {
    result = mLeadTerm;
    return true;
  }

  do {
    if (mQueue.empty())
      return false;
    MultipleWithPos* entry = mQueue.top();
    MATHICGB_ASSERT(entry != 0);

    // remove node from hash table first since we are going to be changing
    // the monomial after this, and if we do that before the hash value will
    // change.
    mHashTable.remove(entry->node);

    // extract information into mLeadTerm
    mLeadTerm.monom.swap(entry->current);
    entry->node->monom = entry->current;
    mLeadTerm.coeff = entry->node->coeff;

    // remove old monomial from hash table and insert next
    MATHICGB_ASSERT(entry->pos != entry->end);
    while (true) {
      ++entry->pos;
      if (entry->pos == entry->end) {
        mQueue.pop();
        entry->destroy(mRing);
        mPool.free(entry);
        break;
      }
      term t;
      t.monom = entry->current;
      entry->computeCurrent(mRing, t.monom);
      entry->currentCoefficient(mRing, t.coeff);

      std::pair<bool, PolyHashTable::node*> p = mHashTable.insert(t);
      if (p.first) {
        entry->node = p.second;
        mQueue.decreaseTop(entry);
        break;
      }
    }
  } while (mRing.coefficientIsZero(mLeadTerm.coeff));

  result = mLeadTerm;
  mLeadTermKnown = true;
  return true;
}

template<template<typename> class Q>
void ReducerHashPack<Q>::removeLeadTerm()
{
  if (!mLeadTermKnown) {
    const_term dummy;
    leadTerm(dummy);
  }
  mLeadTermKnown = false;
}

template<template<typename> class Q>
void ReducerHashPack<Q>::insertEntry(MultipleWithPos* entry) {
  MATHICGB_ASSERT(entry != 0);
  for (; entry->pos != entry->end; ++entry->pos) {
    term t;
    t.monom = entry->current;
    entry->computeCurrent(mRing, t.monom);
    entry->currentCoefficient(mRing, t.coeff);

    std::pair<bool, PolyHashTable::node*> p = mHashTable.insert(t);
    if (p.first) {
      mLeadTermKnown = false;
      entry->node = p.second;
      mQueue.push(entry);
      return;
    }
  }
  entry->destroy(mRing);
  mPool.free(entry);
}

template<template<typename> class Q>
void ReducerHashPack<Q>::resetReducer()
{
  mLeadTermKnown = false;
  MonomialFree freeer(mRing);
#if 0
  //TODO: reinstate these lines, once Geobucket, TourTree and Heap can all handle them
  mQueue.forAll(freeer);
  mQueue.clear();
#endif
  mHashTable.reset();
}

template<template<typename> class Q>
size_t ReducerHashPack<Q>::getMemoryUse() const
{
  return mQueue.getMemoryUse() +
    mPool.getMemoryUse() +
    mHashTable.getMemoryUse();
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
