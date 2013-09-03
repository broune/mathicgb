// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "BjarkeGeobucket.hpp"

MATHICGB_NAMESPACE_BEGIN

BjarkeGeobucket::BjarkeGeobucket(const PolyRing *R0)
  : Reducer(),
    R_(R0),
    H_(R0,10),
    mNodeCount(0),
    G_(Configuration(R0,4,1))
{
}

BjarkeGeobucket::~BjarkeGeobucket()
{
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
void BjarkeGeobucket::insertTail(const_term multiplier, const Poly *g1)
{
  if (g1->nTerms() <= 1) return;

  HashPoly M;
  H_.insert(multiplier, ++(g1->begin()), g1->end(), M);

  if (!M.empty())
    {
      G_.push(M.begin(),M.end());
      mNodeCount += M.size();
    }

  stats_n_inserts++;
  stats_n_compares += G_.getConfiguration().getComparisons();
  G_.getConfiguration().resetComparisons();
}

void BjarkeGeobucket::insert(monomial multiplier, const Poly *g1)
{
  HashPoly M;

  H_.insert(multiplier, g1->begin(), g1->end(), M);

  if (!M.empty())
    {
      G_.push(M.begin(),M.end());
      mNodeCount += M.size();
    }

  stats_n_inserts++;
  stats_n_compares += G_.getConfiguration().getComparisons();
  G_.getConfiguration().resetComparisons();
}

bool BjarkeGeobucket::findLeadTerm(const_term &result)
{
  while (!G_.empty())
    {
      if (H_.popTerm(G_.top(), result.coeff, result.monom))
        // returns true if G_.top() is not the zero element
        return true;
      G_.pop();
      mNodeCount--;
    }
  return false;
}

void BjarkeGeobucket::removeLeadTerm()
// returns true if there is a term to extract
{
  G_.pop();
  mNodeCount--;
}

void BjarkeGeobucket::value(Poly &result)
// keep extracting lead term until done
{
  const_term t;
  while (findLeadTerm(t))
    {
      result.appendTerm(t.coeff, t.monom);
      G_.pop();
      mNodeCount--;
    }
  resetReducer();
}

void BjarkeGeobucket::resetReducer()
{
  const_term t;
  while (findLeadTerm(t))
    {
      G_.pop();
      mNodeCount--;
    }
  H_.reset();
  // how to reset G_ ?
}

size_t BjarkeGeobucket::getMemoryUse() const
{
  size_t result = H_.getMemoryUse();
  result += G_.getMemoryUse();
  //  std::cerr << "[reducer: " << H_.getMemoryUse() << " " << G_.getMemoryUse() << "]" << std::endl;
  return result;
}

MATHICGB_NAMESPACE_END
