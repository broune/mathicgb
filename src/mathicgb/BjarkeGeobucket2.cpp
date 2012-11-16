// Copyright 2011 Michael E. Stillman
#include "stdinc.h"
#include "BjarkeGeobucket2.hpp"

#include <iostream>

BjarkeGeobucket2::BjarkeGeobucket2(const PolyRing *R0):
  mRing(*R0),
  mHashTableOLD(R0, 10),
  mHeap(GeoConfiguration(*R0, 4, 1)),
  mHashTable(BjarkeGeobucket2Configuration(*R0), 10) {
}

void BjarkeGeobucket2::insert(Poly::const_iterator first, 
                              Poly::const_iterator last,
                              std::vector<node*> &result)
{
  for (Poly::const_iterator i = first; i != last; ++i)
    {
      monomial monomspace = mRing.allocMonomial(mArena);
      mRing.monomialCopy(i.getMonomial(), monomspace);
      std::pair<bool, node*> found = mHashTable.insert(monomspace, i.getCoefficient());
      if (found.first)
        {
          // remove the monomial.  It should be at the top of the mArena arena.
          mRing.freeTopMonomial(mArena,monomspace);
          mRing.coefficientAddTo(found.second->value(), i.getCoefficient());
        }
      else
        {
          result.push_back(found.second);
        }
    }
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
void BjarkeGeobucket2::insertTail(const_term multiplier, const Poly *g1)
{
  MATHICGB_ASSERT(g1 != 0);
  MATHICGB_ASSERT(g1->termsAreInDescendingOrder());

  if (g1->nTerms() <= 1)
    return;

  HashPoly M;
  mHashTableOLD.insert(multiplier, ++(g1->begin()), g1->end(), M);

  if (!M.empty())
    mHeap.push(M.begin(),M.end());

  stats_n_inserts++;
  stats_n_compares += mHeap.getConfiguration().getComparisons();
  mHeap.getConfiguration().resetComparisons();
}

void BjarkeGeobucket2::insert(monomial multiplier, const Poly *g1)
{
  MATHICGB_ASSERT(g1 != 0);
  MATHICGB_ASSERT(g1->termsAreInDescendingOrder());

  HashPoly M;
  mHashTableOLD.insert(multiplier, g1->begin(), g1->end(), M);
  if (!M.empty())
    mHeap.push(M.begin(),M.end());

  stats_n_inserts++;
  stats_n_compares += mHeap.getConfiguration().getComparisons();
  mHeap.getConfiguration().resetComparisons();
}

bool BjarkeGeobucket2::leadTerm(const_term &result)
{
  while (!mHeap.empty())
    {
      if (mHashTableOLD.popTerm(mHeap.top(), result.coeff, result.monom))
        // returns true if mHeap.top() is not the zero element
        return true;
      mHeap.pop();
    }
  return false;
}

void BjarkeGeobucket2::removeLeadTerm()
// returns true if there is a term to extract
{
  mHeap.pop();
}

void BjarkeGeobucket2::value(Poly &result)
// keep extracting lead term until done
{
  const_term t;
  while (leadTerm(t))
    {
      result.appendTerm(t.coeff, t.monom);
      mHeap.pop();
    }
  resetReducer();
}

void BjarkeGeobucket2::resetReducer()
{
  const_term t;
  while (leadTerm(t))
    {
      mHeap.pop();
    }
  mHashTableOLD.reset();
  // how to reset mHeap ?
}

size_t BjarkeGeobucket2::getMemoryUse() const
{
  size_t result = TypicalReducer::getMemoryUse();
  result += mHashTableOLD.getMemoryUse();
  result += mHeap.getMemoryUse();
  result += mHashTable.memoryUse();
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
