// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include "BjarkeGeobucket.hpp"

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

  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
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

  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
}

void BjarkeGeobucket::insert(monomial multiplier, const Poly *g1)
{
  HashPoly M;

  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());

  H_.insert(multiplier, g1->begin(), g1->end(), M);

  if (!M.empty())
    {
      G_.push(M.begin(),M.end());
      mNodeCount += M.size();
    }

  stats_n_inserts++;
  stats_n_compares += G_.getConfiguration().getComparisons();
  G_.getConfiguration().resetComparisons();

  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
}

bool BjarkeGeobucket::findLeadTerm(const_term &result)
{
  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
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

  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
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

  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
  resetReducer();
}

void BjarkeGeobucket::resetReducer()
{
  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
  const_term t;
  while (findLeadTerm(t))
    {
      G_.pop();
      mNodeCount--;
    }
  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
  H_.reset();
  MATHICGB_ASSERT(mNodeCount == H_.getNodeCount());
  // how to reset G_ ?
}

size_t BjarkeGeobucket::getMemoryUse() const
{
  size_t result = H_.getMemoryUse();
  result += G_.getMemoryUse();
  //  std::cerr << "[reducer: " << H_.getMemoryUse() << " " << G_.getMemoryUse() << "]" << std::endl;
  return result;
}

void BjarkeGeobucket::dump() const
{
  H_.dump(0);
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
