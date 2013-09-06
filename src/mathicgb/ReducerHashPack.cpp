// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "ReducerHashPack.hpp"

#include "TypicalReducer.hpp"
#include "ReducerHelper.hpp"
#include "PolyHashTable.hpp"
#include <mathic.h>
#include <memtailor.h>

MATHICGB_NAMESPACE_BEGIN

void reducerHashPackDependency() {}

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
    PolyHashTable::Node* node;

    void destroy(const PolyRing& ring);
  };

  class Configuration : public ReducerHelper::PlainConfiguration {
  public:
    typedef MultipleWithPos* Entry;
    Configuration(const PolyRing& ring) : PlainConfiguration(ring) {}
    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monoid().lessThan(a->node->mono(), b->node->mono());
    }
  };

  class MonomialFree;

  void insertEntry(MultipleWithPos* entry);

  const PolyRing& mRing;
  Queue<Configuration> mQueue;
  PolyHashTable mHashTable;
  memt::BufferPool mPool;
};

template<template<typename> class Q>
ReducerHashPack<Q>::ReducerHashPack(const PolyRing& ring):
  mRing(ring),
  mQueue(Configuration(ring)),
  mHashTable(ring),
  mPool(sizeof(MultipleWithPos))
{}

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
  node(0) {}

template<template<typename> class Q>
void ReducerHashPack<Q>::MultipleWithPos::destroy(const PolyRing& ring) {
  ring.freeMonomial(const_cast<ConstMonomial&>(multiple.monom).castAwayConst());

  // Call the destructor to destruct the iterators into std::vector.
  // In debug mode MSVC puts those in a linked list and the destructor
  // has to be called since it takes an iterator off the list. We had
  // memory corruption problems before doing this.
  this->~MultipleWithPos();
}

template<template<typename> class Q>
bool ReducerHashPack<Q>::leadTerm(const_term& result) {
  while (!mQueue.empty()) {
    MultipleWithPos* entry = mQueue.top();
    MATHICGB_ASSERT(entry != 0);

    if (!mRing.coefficientIsZero(entry->node->value())) {
      result.coeff = entry->node->value();
      result.monom = Monoid::toOld(entry->node->mono());
      return true;
    }
    removeLeadTerm();
  }
  return false;
}

template<template<typename> class Q>
void ReducerHashPack<Q>::removeLeadTerm() {
  MATHICGB_ASSERT(!mQueue.empty());

  MultipleWithPos* entry = mQueue.top();
  MATHICGB_ASSERT(entry != 0);

  // remove node from hash table first since we are going to be changing
  // the monomial after this, and if we do that before the hash value will
  // change.
  mHashTable.remove(entry->node);

  MATHICGB_ASSERT(entry->pos != entry->end);
  while (true) {
    ++entry->pos;
    if (entry->pos == entry->end) {
      mQueue.pop();
      entry->destroy(mRing);
      mPool.free(entry);
      break;
    }
   
    const auto p = mHashTable.insertProduct
      (entry->multiple, entry->pos.term());
    if (p.second) {
      entry->node = p.first;
      mQueue.decreaseTop(entry);
      break;
    }
  }
}

template<template<typename> class Q>
void ReducerHashPack<Q>::insertEntry(MultipleWithPos* entry) {
  MATHICGB_ASSERT(entry != 0);
  for (; entry->pos != entry->end; ++entry->pos) {
    const auto p = mHashTable.insertProduct
      (entry->multiple, entry->pos.term());
    if (p.second) {
      entry->node = p.first;
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
  MonomialFree freeer(mRing);
  mQueue.forAll(freeer);
  mQueue.clear();
  mHashTable.clear();
}

template<template<typename> class Q>
size_t ReducerHashPack<Q>::getMemoryUse() const
{
  return mQueue.getMemoryUse() +
    mPool.getMemoryUse() +
    mHashTable.getMemoryUse();
}

Reducer::Registration r6(
  "TourHashPack",
  Reducer::Reducer_TourTree_Hashed_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHashPack<mic::TourTree>>(ring);
  }
);
 
Reducer::Registration r12(
  "HeapHashPack",
  Reducer::Reducer_Heap_Hashed_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHashPack<mic::Heap>>(ring);
  }
);
Reducer::Registration r18(
  "GeoHashPack",
  Reducer::Reducer_Geobucket_Hashed_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHashPack<mic::Geobucket>>(ring);
  }
);

MATHICGB_NAMESPACE_END
