// Copyright 2011 Michael E. Stillman

#include <iostream>
#include "stdinc.h"
#include "BjarkeGeobucket2.hpp"

extern int tracingLevel;

BjarkeGeobucket2::BjarkeGeobucket2(const PolyRing *R0)
  : Reducer(),
    mRing(*R0),
    mHashTableOLD(R0,10),
    mNodeCount(0),
    mHeap(Configuration(*R0,4,1))
{
}

BjarkeGeobucket2::~BjarkeGeobucket2()
{
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
void BjarkeGeobucket2::insertTail(const_term multiplier, const Poly *g1)
{
  if (g1->nTerms() <= 1) return;

  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
  HashPoly M;
  mHashTableOLD.insert(multiplier, ++(g1->begin()), g1->end(), M);

  if (!M.empty())
    {
      mHeap.push(M.begin(),M.end());
      mNodeCount += M.size();
    }

  stats_n_inserts++;
  stats_n_compares += mHeap.getConfiguration().getComparisons();
  mHeap.getConfiguration().resetComparisons();

  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
}

void BjarkeGeobucket2::insert(monomial multiplier, const Poly *g1)
{
  HashPoly M;

  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());

  mHashTableOLD.insert(multiplier, g1->begin(), g1->end(), M);

  if (!M.empty())
    {
      mHeap.push(M.begin(),M.end());
      mNodeCount += M.size();
    }

  stats_n_inserts++;
  stats_n_compares += mHeap.getConfiguration().getComparisons();
  mHeap.getConfiguration().resetComparisons();

  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
}

bool BjarkeGeobucket2::findLeadTerm(const_term &result)
{
  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
  while (!mHeap.empty())
    {
      if (mHashTableOLD.popTerm(mHeap.top(), result.coeff, result.monom))
        // returns true if mHeap.top() is not the zero element
        return true;
      mHeap.pop();
      mNodeCount--;
    }
  return false;
}

void BjarkeGeobucket2::removeLeadTerm()
// returns true if there is a term to extract
{
  mHeap.pop();
  mNodeCount--;

  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
}

void BjarkeGeobucket2::value(Poly &result)
// keep extracting lead term until done
{
  const_term t;
  while (findLeadTerm(t))
    {
      result.appendTerm(t.coeff, t.monom);
      mHeap.pop();
      mNodeCount--;
    }

  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
  resetReducer();
}

void BjarkeGeobucket2::resetReducer()
{
  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
  const_term t;
  while (findLeadTerm(t))
    {
      mHeap.pop();
      mNodeCount--;
    }
  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
  mHashTableOLD.reset();
  ASSERT(mNodeCount == mHashTableOLD.getNodeCount());
  // how to reset mHeap ?
}

size_t BjarkeGeobucket2::getMemoryUse() const
{
  size_t result = mHashTableOLD.getMemoryUse();
  result += mHeap.getMemoryUse();
  //  std::cerr << "[reducer: " << mHashTableOLD.getMemoryUse() << " " << mHeap.getMemoryUse() << "]" << std::endl;
  return result;
}

void BjarkeGeobucket2::dump() const
{
  mHashTableOLD.dump(0);
}

// Local Variables:
// compile-command: "make -C $MATHIC/mathicgb "
// indent-tabs-mode: nil
// End:
