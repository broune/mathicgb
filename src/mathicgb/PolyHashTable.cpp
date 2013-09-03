// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "PolyHashTable.hpp"

#include <mathic.h>
#include <iostream>
#include <cmath>

MATHICGB_NAMESPACE_BEGIN

PolyHashTable::PolyHashTable(const PolyRing *R, int nbits)
  : mRing(*R),
    mHashMask((static_cast<size_t>(1) << nbits)-1),
    mTableSize(static_cast<size_t>(1) << nbits),
    mLogTableSize(nbits),
    mNodeCount(0),
    mBinCount(0),
    mMaxCountBeforeRebuild(0)
{
  mHashTable.resize(mTableSize);

  mMonomialSize = R->maxMonomialSize() * sizeof(exponent);
  // set each entry of mHashTable to null

  reset();
}

void PolyHashTable::reset()
{
  // The following no longer needs to be done

  // Clear the table, and memory areas.

#if 0
  MATHICGB_ASSERT(mNodeCount != 0);
  for (size_t count = 0; count < mTableSize; ++count)
    {
      if ( mHashTable[count] != 0)
        std::cerr << "error: hash table is not zero on reset" << std::endl;
    }
#endif

  mArena.freeAllAllocs();

  mBinCount = 0;
  mNodeCount = 0;
}

void PolyHashTable::resize(size_t new_nbits)
// Don't change the nodes, table, but do recreate mHashTable
{
  // Make a new vector of node *'s.
  // swap the two.
  // Loop through each one, reinserting the node into the proper bin.

  //  std::cout << "resizing PolyHashTable to " << new_nbits << " bits" << " count=" << mNodeCount << std::endl;
  size_t const old_table_size = mTableSize;
  mTableSize = static_cast<size_t>(1) << new_nbits;
  mLogTableSize = new_nbits;
  mHashMask = mTableSize-1;
  MonomialArray old_table(mTableSize);

  std::swap(old_table, mHashTable);
  mBinCount = 0;
  
  for (size_t i = 0; i < old_table_size; ++i)
    {
      node *p = old_table[i];
      while (p != 0)
        {
          node *q = p;
          p = p->next();
          q->next() = 0;
          // Reinsert node.  We know that it is unique
          const_monomial m = q->monom;
          size_t hashval = mRing.monomialHashValue(m) & mHashMask;
          node *r = mHashTable[hashval];
          if (r == 0) 
            {
              mBinCount++;
              q->next() = r;
              mHashTable[hashval] = q;
            }
          else
            {
              // put it at the end
              for ( ; r->next() != 0; r = r->next()) { }
              r->next() = q;
            }
        }
    }

  // todo: consider if this can overflow or something else nasty might happen
  const double threshold = 0.1;
  mMaxCountBeforeRebuild =
    static_cast<size_t>(std::floor(mTableSize * threshold));
}

PolyHashTable::node * PolyHashTable::makeNode(coefficient coeff, const_monomial monom)
{
  mNodeCount++;
  node *q = static_cast<node *>(mArena.allocObjectNoCon<node>());
  q->next() = 0;
  q->monom = monom; 
  mRing.coefficientSet(q->coeff, coeff);
  return q;
}

bool PolyHashTable::lookup_and_insert(const_monomial m, coefficient val, node *& result)
// Returns true if m is in the table, else inserts m into the hash table (as is, without copying it)
{
  size_t fullHashVal = mRing.monomialHashValue(m);
  size_t hashval = fullHashVal & mHashMask;

  MATHICGB_ASSERT(hashval < mHashTable.size());
  node *tmpNode = mHashTable[hashval];
  if (tmpNode == 0) {
    result = makeNode(val, m);
    mHashTable[hashval] = result;
  } else {
    while (true) {
      if (mRing.monomialHashValue(tmpNode->monom) == fullHashVal && mRing.monomialEQ(m, tmpNode->monom)) {
        mRing.coefficientAddTo(tmpNode->coeff, val);
        result = tmpNode;
        return true;
      }
      if (tmpNode->next() == 0) {
        result = makeNode(val, m);
        tmpNode->next() = result;
        break;
      }
      tmpNode = tmpNode->next();
    }
  }

  if (mNodeCount > mMaxCountBeforeRebuild)
    resize(mLogTableSize + 2);  // increase by a factor of 4??

  return false;
}

void PolyHashTable::insert(
  Poly::const_iterator first, 
  Poly::const_iterator last,
  MonomialArray &result
) {
  for (auto i = first; i != last; ++i) {
    monomial monomspace = mRing.allocMonomial(mArena);
    node* p;
    mRing.monomialCopy(i.getMonomial(), monomspace);
    bool found = lookup_and_insert(monomspace, i.getCoefficient(), p);
    if (found)
      mRing.freeTopMonomial(mArena,monomspace);
    else
      result.push_back(p);
  }
}

void PolyHashTable::insert(
  const_term multiplier, 
  Poly::const_iterator first, 
  Poly::const_iterator last,
  MonomialArray &result
) {
  for (auto i = first; i != last; ++i) {
    monomial monomspace = mRing.allocMonomial(mArena);
    coefficient c;
    mRing.coefficientSet(c, multiplier.coeff);
    node* p;
    mRing.monomialMult(multiplier.monom, i.getMonomial(), monomspace);
    mRing.coefficientMultTo(c, i.getCoefficient());
    bool found = lookup_and_insert(monomspace, c, p);
    if (found)
      mRing.freeTopMonomial(mArena,monomspace);
    else
      result.push_back(p);
  }
}

void PolyHashTable::insert(
  const_monomial multiplier, 
  Poly::const_iterator first, 
  Poly::const_iterator last,
  MonomialArray& result
) {
  for (Poly::const_iterator i = first; i != last; ++i) {
    monomial monomspace = mRing.allocMonomial(mArena);
    node* p;
    mRing.monomialMult(multiplier, i.getMonomial(), monomspace);
    bool found = lookup_and_insert(monomspace, i.getCoefficient(), p);
    if (found)
      mRing.freeTopMonomial(mArena,monomspace);
    else
      result.push_back(p);
  }
}

std::pair<bool, PolyHashTable::node*>
PolyHashTable::insert(const_term termToInsert) {
  node* n;
  bool alreadyInThere = lookup_and_insert
    (termToInsert.monom, termToInsert.coeff, n);
  return std::make_pair(!alreadyInThere, n);
}

void PolyHashTable::unlink(node* p)
{
  mNodeCount--;
  size_t const hashval = mRing.monomialHashValue(p->monom) & mHashMask;

  node head;
  node* tmpNode = mHashTable[hashval];
  head.next() = tmpNode;
  for (node* q = &head; q->next() != 0; q = q->next()) {
    if (q->next() == p) {
      q->next() = p->next();
      mHashTable[hashval] = head.next();
      return;
    }
  }
  // If we get here, then the node is not at its supposed hash value.
  // That probably means either that the node has been deleted twice
  // or that the value in the node changed so that its hash value
  // changed. That is not allowed.
  MATHICGB_ASSERT(false);
}

void PolyHashTable::remove(node* n) {
  unlink(n);
}

bool PolyHashTable::popTerm(node *p, coefficient &result_coeff, const_monomial &result_monom)
{
  unlink(p);
  if (!mRing.coefficientIsZero(p->coeff))
    {
      result_coeff = p->coeff;
      result_monom = p->monom;
      return true;
    }
  return false;
}

size_t PolyHashTable::getMemoryUse() const
{
  size_t result = mHashTable.capacity() * sizeof(node *);
  result += mArena.getMemoryUse();
  return result;
}

MATHICGB_NAMESPACE_END
