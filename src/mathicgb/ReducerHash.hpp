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
  ~ReducerHash();

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
    typedef PolyHashTable::node * Entry;

    Configuration(const PolyRing& ring): 
      PlainConfiguration(ring),
      mComparisonCount(0) {}

    CompareResult compare(const Entry& a, const Entry& b) const {
      ++mComparisonCount;
      return ring().monomialLT(a->monom, b->monom);
    }

    unsigned long long getComparisonCount() const {return mComparisonCount;}

    void resetComparisonCount() const {mComparisonCount = 0;}
    
  private:
    mutable unsigned long long mComparisonCount;
  };
  
private:
  const PolyRing &mRing;
  PolyHashTable mHashTable;
  Queue<Configuration> mQueue;

  // Number of (distinct) monomials in mQueue.  
  // Statistics and debugging use only
  size_t mNodeCount;  
};

template<template<typename> class Q>
ReducerHash<Q>::ReducerHash(const PolyRing &ring):
  mRing(ring),
  mHashTable(&ring,10),
  mQueue(Configuration(ring)),
  mNodeCount(0) {
}

template<template<typename> class Q>
ReducerHash<Q>::~ReducerHash()
{
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////

template<template<typename> class Q>
void ReducerHash<Q>::insertTail(const_term multiplier, const Poly *g1)
{
  if (g1->nTerms() <= 1) return;

  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
  PolyHashTable::MonomialArray M;
  mHashTable.insert(multiplier, ++(g1->begin()), g1->end(), M);

  if (!M.empty()) {
    mQueue.push(M.begin(),M.end());
    mNodeCount += M.size();
  }

  stats_n_inserts++;
  stats_n_compares += mQueue.getConfiguration().getComparisonCount();
  mQueue.getConfiguration().resetComparisonCount();

  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
}

template<template<typename> class Q>
void ReducerHash<Q>::insert(monomial multiplier, const Poly *g1)
{
  PolyHashTable::MonomialArray M;

  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());

  mHashTable.insert(multiplier, g1->begin(), g1->end(), M);

  if (!M.empty())
    {
      mQueue.push(M.begin(),M.end());
#if 0
      for (PolyHashTable::MonomialArray::const_iterator a = M.begin(); a != M.end(); ++a)
        mQueue.push(*a);
#endif
      mNodeCount += M.size();
    }

  stats_n_inserts++;
  stats_n_compares += mQueue.getConfiguration().getComparisonCount();
  mQueue.getConfiguration().resetComparisonCount();

  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
}

template<template<typename> class Q>
bool ReducerHash<Q>::leadTerm(const_term &result)
{
  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
  while (!mQueue.empty())
    {
      if (mHashTable.popTerm(mQueue.top(), result.coeff, result.monom))
        // returns true if mQueue.top() is not the zero element
        return true;
      mQueue.pop();
      mNodeCount--;
    }
  return false;
}

template<template<typename> class Q>
void ReducerHash<Q>::removeLeadTerm()
// returns true if there is a term to extract
{
  mQueue.pop();
  mNodeCount--;

  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
}

template<template<typename> class Q>
void ReducerHash<Q>::resetReducer()
{
  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
  const_term t;
  while (leadTerm(t))
    {
      mQueue.pop();
      mNodeCount--;
    }
  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
  mHashTable.reset();
  MATHICGB_ASSERT(mNodeCount == mHashTable.getNodeCount());
  // how to reset mQueue ?
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
