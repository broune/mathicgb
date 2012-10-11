// Copyright 2011 Michael E. Stillman
#include "stdinc.h"
#include "PolyHashTable.hpp"

#include <mathic.h>
#include <iostream>
#include <cmath>


const double PolyHashTable::threshold = 0.1;
const bool AlwaysInsertAtEnd = true;

PolyHashTable::PolyHashTable(const PolyRing *R, int nbits)
  : mRing(*R),
    mHashMask((1 << nbits)-1),
    mTableSize(1 << nbits),
    mLogTableSize(nbits),
    mNodeCount(0),
    mBinCount(0),
    mMaxCountBeforeRebuild(0)
{
  mHashTable.resize(mTableSize);

  mMonomialSize = R->maxMonomialSize() * sizeof(exponent);
  // set each entry of mHashTable to null

  mStats.max_table_size = 0;
  mStats.max_chain_len_ever = 0;
  mStats.n_resets = 0;
  mStats.max_n_nonempty_bins = 0;
  mStats.n_inserts = 0;
  mStats.n_moneq_true = 0;
  mStats.n_moneq_false = 0;
  mStats.n_easy_inserts = 0;

  mStats.n_nodes = 0;
  mStats.n_nonempty_bins = 0;
  reset();
}

void PolyHashTable::reset()
{
  // The following no longer needs to be done

  // Clear the table, and memory areas.

  mStats.n_resets++;

#if 0
  ASSERT(mNodeCount != 0);
  for (size_t count = 0; count < mTableSize; ++count)
    {
      if ( mHashTable[count] != 0)
        std::cerr << "error: hash table is not zero on reset" << std::endl;
    }
#endif

  mArena.freeAllAllocs();

  mBinCount = 0;
  mNodeCount = 0;
  
  resetStats();
  MATHICGB_SLOW_ASSERT(computeNodeCount() == 0);
}

size_t PolyHashTable::computeNodeCount() const
{
  size_t result = 0;
  for (size_t i=0; i<mTableSize; i++)
    {
      for (node *p=mHashTable[i]; p != 0; p=p->next) result++;
    }
  return result;
}
void PolyHashTable::resize(size_t new_nbits)
// Don't change the nodes, table, but do recreate mHashTable
{
  // Make a new vector of node *'s.
  // swap the two.
  // Loop through each one, reinserting the node into the proper bin.

  //  std::cout << "resizing PolyHashTable to " << new_nbits << " bits" << " count=" << mNodeCount << std::endl;
  ASSERT(computeNodeCount() == mNodeCount);
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
          p = p->next;
          q->next = 0;
          // Reinsert node.  We know that it is unique
          const_monomial m = q->monom;
          size_t hashval = mRing.monomialHashValue(m) & mHashMask;
          node *r = mHashTable[hashval];
          if (r == 0 || !AlwaysInsertAtEnd) 
            {
              mBinCount++;
              q->next = r;
              mHashTable[hashval] = q;
            }
          else
            {
              // put it at the end
              for ( ; r->next != 0; r = r->next) { }
              r->next = q;
            }
        }
    }

  mStats.max_table_size = mTableSize;

  // todo: consider if this can overflow or something else nasty might happen
  mMaxCountBeforeRebuild =
    static_cast<size_t>(std::floor(mTableSize * threshold));

  ASSERT(computeNodeCount() == mNodeCount);
}

PolyHashTable::node * PolyHashTable::makeNode(coefficient coeff, const_monomial monom)
{
  mNodeCount++;
  if (mNodeCount > mStats.n_nodes)
    mStats.n_nodes = mNodeCount;
  if (mBinCount > mStats.n_nonempty_bins)
    mStats.n_nonempty_bins = mBinCount;
  node *q = static_cast<node *>(mArena.allocObjectNoCon<node>());
  q->next = 0;
  q->monom = monom; 
  mRing.coefficientSet(q->coeff, coeff);
  return q;
}

bool PolyHashTable::lookup_and_insert(const_monomial m, coefficient val, node *& result)
// Returns true if m is in the table, else inserts m into the hash table (as is, without copying it)
{
  mStats.n_inserts++;

  size_t fullHashVal = mRing.monomialHashValue(m);
  size_t hashval = fullHashVal & mHashMask;

  ASSERT(hashval < mHashTable.size());
  node *tmpNode = mHashTable[hashval];
  if (tmpNode == 0)
    {
      mStats.n_easy_inserts++;
      result = makeNode(val, m);
      mHashTable[hashval] = result;
    }
  else
    {
      // loop through to see if we have it
      size_t chainLength = 0;
      while (true)
        {
          if (mRing.monomialHashValue(tmpNode->monom) == fullHashVal && mRing.monomialEQ(m, tmpNode->monom))
            {
              mStats.n_moneq_true++;
              mRing.coefficientAddTo(tmpNode->coeff, val);
              result = tmpNode;
              return true;
            }
          mStats.n_moneq_false++;
          if (tmpNode->next == 0)
            {
              // time to insert the monomial
              result = makeNode(val, m);
              chainLength++;
              if (AlwaysInsertAtEnd)
                {
                  tmpNode->next = result;
                }
              else
                {
                  result->next = mHashTable[hashval];
                  mHashTable[hashval] = result;
                }
              break;
            }
          tmpNode = tmpNode->next;
          chainLength++;
        }
      if (chainLength > mStats.max_chain_len_ever)
        mStats.max_chain_len_ever = chainLength;
    }


  if (mNodeCount > mMaxCountBeforeRebuild)
    resize(mLogTableSize + 2);  // increase by a factor of 4??

  MATHICGB_SLOW_ASSERT(computeNodeCount() == mNodeCount);

  return false;
}

void PolyHashTable::insert(Poly::const_iterator first, 
                           Poly::const_iterator last,
                           MonomialArray &result)
{
  for (Poly::const_iterator i = first; i != last; ++i)
    {
      monomial monomspace = mRing.allocMonomial(mArena);
      node *p;
      mRing.monomialCopy(i.getMonomial(), monomspace);
      bool found = lookup_and_insert(monomspace, i.getCoefficient(), p);
      if (found)
        {
          // remove the monomial.  It should be at the top of the mArena arena.
          mRing.freeTopMonomial(mArena,monomspace);
        }
      else
        {
          result.push_back(p);
        }
    }
}

void PolyHashTable::insert(const_term multiplier, 
                           Poly::const_iterator first, 
                           Poly::const_iterator last,
                           MonomialArray &result)
{
  for (Poly::const_iterator i = first; i != last; ++i)
    {
      monomial monomspace = mRing.allocMonomial(mArena);
      coefficient c;
      mRing.coefficientSet(c, multiplier.coeff);
      node *p;
      mRing.monomialMult(multiplier.monom, i.getMonomial(), monomspace);
      mRing.coefficientMultTo(c, i.getCoefficient());
      bool found = lookup_and_insert(monomspace, c, p);
      if (found)
        {
          // remove the monomial.  It should be at the top of the mArena arena.
          mRing.freeTopMonomial(mArena,monomspace);
        }
      else
        {
          result.push_back(p);
        }
    }
}

void PolyHashTable::insert(const_monomial multiplier, 
                           Poly::const_iterator first, 
                           Poly::const_iterator last,
                           MonomialArray &result)
{
  for (Poly::const_iterator i = first; i != last; ++i)
    {
      monomial monomspace = mRing.allocMonomial(mArena);
      node *p;
      mRing.monomialMult(multiplier, i.getMonomial(), monomspace);
      bool found = lookup_and_insert(monomspace, i.getCoefficient(), p);
      if (found)
        {
          // remove the monomial. It should be at the top of the mArena arena.
          mRing.freeTopMonomial(mArena,monomspace);
        }
      else
        {
          result.push_back(p);
        }
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
  head.next = tmpNode;
  for (node* q = &head; q->next != 0; q = q->next) {
    if (q->next == p) {
      q->next = p->next;
      mHashTable[hashval] = head.next;
      return;
    }
  }
  // If we get here, then the node is not at its supposed hash value.
  // That probably means either that the node has been deleted twice
  // or that the value in the node changed so that its hash value
  // changed. That is not allowed.
  ASSERT(false);
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

void PolyHashTable::toPoly(const MonomialArray::const_iterator &fbegin,
                           const MonomialArray::const_iterator &fend,
                           Poly &result)
{
  coefficient coeff;
  const_monomial monom;
  for (MonomialArray::const_iterator i = fbegin; i != fend; ++i)
    if (popTerm(*i, coeff, monom))
      result.appendTerm(coeff, monom);
}

void PolyHashTable::toPoly(const MonomialArray &f, Poly &result)
{
  // Here we take the monomials in f.    Find corresponding coeff, and append to result.
  // ASSUMPTION: The monomials of f are in order, AND each appears in the hash table
  toPoly(f.begin(), f.end(), result);
}

void PolyHashTable::fromPoly(const Poly &f, MonomialArray &result)
{
  insert(f.begin(), f.end(), result);
}

size_t PolyHashTable::getMemoryUse() const
{
  size_t result = mHashTable.capacity() * sizeof(node *);
  result += mArena.getMemoryUse();
  return result;
}

void PolyHashTable::resetStats() const
{
  //  mStats.max_chain_len = 0;
  mStats.n_nonempty_bins = 0;
}

void PolyHashTable::getStats(Stats &stats) const
{
  // First we set the values in mStats

  //  mStats.max_chain_len = 0;

#if 0
  mStats.n_nonempty_bins = 0;
  mStats.n_nodes = 0;
  for (size_t i = 0; i<mTableSize; i++)
    {
      if (mHashTable[i] == 0) continue;
      mStats.n_nonempty_bins++;
      size_t chain_len = 0;
      for (node *p = mHashTable[i]; p != 0; p = p->next)
        chain_len++;
      mStats.n_nodes += chain_len;
      if (chain_len > mStats.max_chain_len)
        mStats.max_chain_len = chain_len;
    }

  if (mStats.max_chain_len > mStats.max_chain_len_ever)
    mStats.max_chain_len_ever = mStats.max_chain_len;
#endif

  if (&stats != &mStats)
    stats = mStats;
}

void PolyHashTable::dump(int level) const
{
  mic::ColumnPrinter pr;
  pr.addColumn();
  pr.addColumn(false);
  pr.addColumn(false);

  std::ostream& name = pr[0];
  std::ostream& value = pr[1];
  std::ostream& extra = pr[2];

  name << "PolyHashTable stats:" << '\n';
  value << "\n";
  extra << "\n";

  mRing.displayHashValues();

  name << "# resets:\n";
  value << mic::ColumnPrinter::commafy(mStats.n_resets) << '\n';
  extra << '\n';

  name << "# bins:\n";
  value << mic::ColumnPrinter::commafy(mTableSize) << '\n';
  extra << '\n';

  name << "max # monomials in table:\n";
  value << mic::ColumnPrinter::commafy(mStats.n_nodes) << '\n';
  extra << '\n';

  name << "max # nonempty bins:\n";
  value << mic::ColumnPrinter::commafy(mStats.n_nonempty_bins) << '\n';
  extra << mic::ColumnPrinter::percent(mStats.n_nonempty_bins, mTableSize) << " used\n";

  name << "max chain length ever:\n";
  value << mic::ColumnPrinter::commafy(mStats.max_chain_len_ever) << '\n';
  extra << '\n';

  //  name << "max chain length this computation:\n";
  //  value << mic::ColumnPrinter::commafy(mStats.max_chain_len) << '\n';
  //  extra << '\n';

  name << "# calls to lookup_and_insert:\n";
  value << mic::ColumnPrinter::commafy(mStats.n_inserts) << '\n';
  extra << '\n';

  name << "# easy inserts:\n";
  value << mic::ColumnPrinter::commafy(mStats.n_easy_inserts) << '\n';
  extra << '\n';

  name << "# monomialEQ true calls:\n";
  value << mic::ColumnPrinter::commafy(mStats.n_moneq_true) << '\n';
  extra << '\n';

  name << "# monomialEQ false calls:\n";
  value << mic::ColumnPrinter::commafy(mStats.n_moneq_false) << '\n';
  extra << "(Also number of monomials inserted in populated bins)\n";

  std::cout << pr << std::flush;

  if (level == 0) return;

  for (size_t i = 0; i<mTableSize; i++)
    {
      if (mHashTable[i] == 0) continue;
      std::cout << "bin " << i << ": ";
      Poly f(mRing);
      for (node *p = mHashTable[i]; p != 0; p = p->next)
        f.appendTerm(p->coeff, p->monom);
      f.display(std::cout);
      std::cout << std::endl;
    }

}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
