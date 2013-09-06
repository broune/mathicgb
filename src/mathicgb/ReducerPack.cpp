// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "ReducerPack.hpp"

#include "TypicalReducer.hpp"
#include "ReducerHelper.hpp"
#include <memtailor.h>
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

void reducerPackDependency() {}

template<template<typename> class Queue>
class ReducerPack : public TypicalReducer {
public:
  ReducerPack(const PolyRing& ring):
    mRing(ring),
    mLeadTerm(0, mRing.allocMonomial()),
    mLeadTermKnown(false),
    mQueue(Configuration(ring)),
    mPool(sizeof(MultipleWithPos))
  {}

  virtual ~ReducerPack() {
    resetReducer();
    mRing.freeMonomial(mLeadTerm.monom);
  }

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
  class MonomialFree {
  public:
    MonomialFree(const PolyRing& ring): mRing(ring) {}

    bool proceed(MultipleWithPos* entry) {
      entry->destroy(mRing);
      return true;
    }

  private:
    const PolyRing& mRing;
  };
  
  const PolyRing& mRing;
  term mLeadTerm;
  bool mLeadTermKnown;
  Queue<Configuration> mQueue;
  memt::BufferPool mPool;
};

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

MATHICGB_REGISTER_REDUCER(
  "TourNoDedupPack",
  Reducer_TourTree_NoDedup_Packed,
  make_unique<ReducerPack<mic::TourTree>>(ring)
);

MATHICGB_REGISTER_REDUCER(
  "HeapNoDedupPack",
  Reducer_Heap_NoDedup_Packed,
  make_unique<ReducerPack<mic::Heap>>(ring)
);

MATHICGB_REGISTER_REDUCER(
  "GeoNoDedupPack",
  Reducer_Geobucket_NoDedup_Packed,
  make_unique<ReducerPack<mic::Geobucket>>(ring)
);

MATHICGB_NAMESPACE_END
