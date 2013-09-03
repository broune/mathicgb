// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_REDUCER_HASH_GUARD
#define MATHICGB_REDUCER_HASH_GUARD

#include "TypicalReducer.hpp"
#include "ReducerHelper.hpp"
#include "PolyHashTable.hpp"
#include <memtailor.h>
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

template<template<typename ConfigType> class Queue> class ReducerHash;

template<template<typename> class Queue>
class ReducerHash : public TypicalReducer {
public:
  ReducerHash(const PolyRing &ring);

  virtual std::string description() const { 
    return mQueue.getName() + "-hashed";
  }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term& result);
  void removeLeadTerm();

  size_t getMemoryUse() const;

protected:
  void resetReducer();

public:
  class Configuration : public ReducerHelper::PlainConfiguration {
  public:
    typedef PolyHashTable::Node* Entry;

    Configuration(const PolyRing& ring): PlainConfiguration(ring) {}

    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monoid().lessThan(a->mono(), b->mono());
    }
  };
  
private:
  mutable std::vector<PolyHashTable::Node*> mNodesTmp;
  const PolyRing &mRing;
  PolyHashTable mHashTable;
  Queue<Configuration> mQueue;
};

template<template<typename> class Q>
ReducerHash<Q>::ReducerHash(const PolyRing &ring):
  mRing(ring),
  mHashTable(ring),
  mQueue(Configuration(ring))
{}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////

template<template<typename> class Q>
void ReducerHash<Q>::insertTail(const_term multiplier, const Poly *g1)
{
  if (g1->nTerms() <= 1) return;

  mNodesTmp.clear();
  auto it = g1->begin();
  const auto end = g1->end();
  for (++it; it != end; ++it) {
    auto p = mHashTable.insertProduct(it.term(), multiplier);
    if (p.second)
      mNodesTmp.emplace_back(p.first);
  }
  if (!mNodesTmp.empty())
    mQueue.push(mNodesTmp.begin(), mNodesTmp.end());
}

template<template<typename> class Q>
void ReducerHash<Q>::insert(monomial multiplier, const Poly *g1)
{
  mNodesTmp.clear();
  const auto end = g1->end();
  for (auto it = g1->begin(); it != end; ++it) {
    auto p = mHashTable.insertProduct
      (it.getMonomial(), multiplier, it.getCoefficient());
    if (p.second)
      mNodesTmp.emplace_back(p.first);
  }
  if (!mNodesTmp.empty())
    mQueue.push(mNodesTmp.begin(), mNodesTmp.end());
}

template<template<typename> class Q>
bool ReducerHash<Q>::leadTerm(const_term& result)
{
  while (!mQueue.empty()) {
    const auto top = mQueue.top();
    if (!mRing.coefficientIsZero(top->value())) {
      result.coeff = top->value();
      result.monom = Monoid::toOld(top->mono());
      return true;
    }
    mQueue.pop();
    mHashTable.remove(top);
  }
  return false;
}

template<template<typename> class Q>
void ReducerHash<Q>::removeLeadTerm()
// returns true if there is a term to extract
{
  const auto top = mQueue.top();
  mQueue.pop();
  mHashTable.remove(top);
}

template<template<typename> class Q>
void ReducerHash<Q>::resetReducer()
{
  while (!mQueue.empty()) {
    const auto top = mQueue.top();
    mQueue.pop();
    mHashTable.remove(top);
  }
  mHashTable.clear();
}

template<template<typename> class Q>
size_t ReducerHash<Q>::getMemoryUse() const
{
  size_t result = TypicalReducer::getMemoryUse();
  result += mHashTable.getMemoryUse();
  result += mQueue.getMemoryUse();
  return result;
}

MATHICGB_NAMESPACE_END
#endif
