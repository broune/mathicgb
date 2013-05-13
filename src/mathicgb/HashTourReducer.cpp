// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "HashTourReducer.hpp"

#include <utility>

MATHICGB_NAMESPACE_BEGIN

HashTourReducer::HashTourReducer(const PolyRing& ring):
  mRing(ring),
  mLeadTerm(0, mRing.allocMonomial()),
  mLeadTermKnown(false),
  mQueue(Configuration(ring)),
  mHashTable(&ring, 10),
  mPool(sizeof(MultipleWithPos))
{
}

class HashTourReducer::MonomialFree
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

HashTourReducer::~HashTourReducer()
{
  resetReducer();
  mRing.freeMonomial(mLeadTerm.monom);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
void HashTourReducer::insertTail(const_term multiple, const Poly* poly)
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

void HashTourReducer::insert(monomial multiple, const Poly* poly)
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

HashTourReducer::MultipleWithPos::MultipleWithPos
(const Poly& poly, const_term multiple):
  pos(poly.begin()),
  end(poly.end()),
  multiple(allocTerm(poly.ring(), multiple)),
  current(poly.ring().allocMonomial()),
  node(0) {}

void HashTourReducer::MultipleWithPos::
computeCurrent(const PolyRing& ring, monomial current) {
  ring.monomialMult(multiple.monom, pos.getMonomial(), current);  
}

void HashTourReducer::MultipleWithPos::currentCoefficient
(const PolyRing& ring, coefficient& coeff) {
  ring.coefficientMult(multiple.coeff, pos.getCoefficient(), coeff);
}

void HashTourReducer::MultipleWithPos::destroy(const PolyRing& ring) {
  ring.freeMonomial(current);
  ring.freeMonomial(const_cast<ConstMonomial&>(multiple.monom).castAwayConst());

  // Call the destructor to destruct the iterators into std::vector.
  // In debug mode MSVC puts those in a linked list and the destructor
  // has to be called since it takes an iterator off the list. We had
  // memory corruption problems before doing this.
  this->~MultipleWithPos();
}

bool HashTourReducer::leadTerm(const_term& result)
{
  if (mLeadTermKnown) {
    result = mLeadTerm;
    return true;
  }

  do {
    if (mQueue.empty())
      return false;
    MultipleWithPos* entry = mQueue.top();
    MATHICGB_ASSERT(entry != 0);

    // remove node from hash table first since we are going to be changing
    // the monomial after this, and if we do that before the hash value will
    // change.
    mHashTable.remove(entry->node);

    // extract information into mLeadTerm
    mLeadTerm.monom.swap(entry->current);
    entry->node->monom = entry->current;
    mLeadTerm.coeff = entry->node->coeff;

    // remove old monomial from hash table and insert next
    MATHICGB_ASSERT(entry->pos != entry->end);
    while (true) {
      ++entry->pos;
      if (entry->pos == entry->end) {
        mQueue.pop();
        entry->destroy(mRing);
        mPool.free(entry);
        break;
      }
      term t;
      t.monom = entry->current;
      entry->computeCurrent(mRing, t.monom);
      entry->currentCoefficient(mRing, t.coeff);

      std::pair<bool, PolyHashTable::node*> p = mHashTable.insert(t);
      if (p.first) {
        entry->node = p.second;
        mQueue.decreaseTop(entry);
        break;
      }
    }
  } while (mRing.coefficientIsZero(mLeadTerm.coeff));

  result = mLeadTerm;
  mLeadTermKnown = true;
  return true;
}

void HashTourReducer::removeLeadTerm()
{
  if (!mLeadTermKnown) {
    const_term dummy;
    leadTerm(dummy);
  }
  mLeadTermKnown = false;
}

void HashTourReducer::insertEntry(MultipleWithPos* entry) {
  MATHICGB_ASSERT(entry != 0);
  for (; entry->pos != entry->end; ++entry->pos) {
    term t;
    t.monom = entry->current;
    entry->computeCurrent(mRing, t.monom);
    entry->currentCoefficient(mRing, t.coeff);

    std::pair<bool, PolyHashTable::node*> p = mHashTable.insert(t);
    if (p.first) {
      mLeadTermKnown = false;
      entry->node = p.second;
      mQueue.push(entry);
      return;
    }
  }
  entry->destroy(mRing);
  mPool.free(entry);
}

void HashTourReducer::value(Poly &result)
{
  const_term t;
  while (leadTerm(t)) {
    result.appendTerm(t.coeff, t.monom);
    removeLeadTerm();
  }
  resetReducer();
}

void HashTourReducer::resetReducer()
{
  mLeadTermKnown = false;
  MonomialFree freeer(mRing);
  mQueue.forAll(freeer);
  mQueue.clear();
  mHashTable.reset();
}

size_t HashTourReducer::getMemoryUse() const
{
  return
    TypicalReducer::getMemoryUse() +
    mQueue.getMemoryUse() +
    mPool.getMemoryUse() +
    mHashTable.getMemoryUse();
}

void HashTourReducer::dump() const
{
}

MATHICGB_NAMESPACE_END
