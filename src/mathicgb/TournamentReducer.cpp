// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "TournamentReducer.hpp"

#include <utility>

MATHICGB_NAMESPACE_BEGIN

TournamentReducer::TournamentReducer(const PolyRing& ring):
  mRing(ring),
  mLeadTerm(0, mRing.allocMonomial()),
  mLeadTermKnown(false),
  mQueue(Configuration(ring)),
  mPool(sizeof(MultipleWithPos))
{
}

class TournamentReducer::MonomialFree
{
public:
  MonomialFree(const PolyRing& ring): mRing(ring) {}

  bool proceed(MultipleWithPos* entry)
  {
    entry->destroy(mRing);
    return true;
  }
private:
  const PolyRing& mRing;
};

TournamentReducer::~TournamentReducer()
{
  resetReducer();
  mRing.freeMonomial(mLeadTerm.monom);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////
void TournamentReducer::insertTail(const_term multiple, const Poly* poly)
{
  if (poly->nTerms() <= 1)
    return;
  mLeadTermKnown = false;

  MultipleWithPos* entry = new (mPool.alloc()) MultipleWithPos(*poly, multiple);
  ++entry->pos;
  entry->computeCurrent(poly->ring());
  mQueue.push(entry);
}

void TournamentReducer::insert(monomial multiple, const Poly* poly)
{
  if (poly->isZero())
    return;
  mLeadTermKnown = false;

  term termMultiple(1, multiple);
  MultipleWithPos* entry = new (mPool.alloc()) MultipleWithPos(*poly, termMultiple);
  entry->computeCurrent(poly->ring());
  mQueue.push(entry);
}

namespace {
  const_term allocTerm(const PolyRing& ring, const_term term) {
    monomial mono = ring.allocMonomial();
    ring.monomialCopy(term.monom, mono);
    return const_term(term.coeff, mono);
  }
}

TournamentReducer::MultipleWithPos::MultipleWithPos
(const Poly& poly, const_term multiple):
  pos(poly.begin()),
  end(poly.end()),
  multiple(allocTerm(poly.ring(), multiple)),
  current(poly.ring().allocMonomial()) {}

void TournamentReducer::MultipleWithPos::computeCurrent(const PolyRing& ring) {
  ring.monomialMult(multiple.monom, pos.getMonomial(), current);  
}

void TournamentReducer::MultipleWithPos::currentCoefficient
(const PolyRing& ring, coefficient& coeff) {
  ring.coefficientMult(multiple.coeff, pos.getCoefficient(), coeff);
}

void TournamentReducer::MultipleWithPos::destroy(const PolyRing& ring) {
  ring.freeMonomial(current);
  ring.freeMonomial(const_cast<ConstMonomial&>(multiple.monom).castAwayConst());

  // Call the destructor to destruct the iterators into std::vector.
  // In debug mode MSVC puts those in a linked list and the destructor
  // has to be called since it takes an iterator off the list. We had
  // memory corruption problems before doing this.
  this->~MultipleWithPos();
}

bool TournamentReducer::leadTerm(const_term& result)
{
  if (mLeadTermKnown) {
    result = mLeadTerm;
    return true;
  }

  do {
    if (mQueue.empty())
      return false;
    MultipleWithPos* entry = mQueue.top();
    mLeadTerm.monom.swap(entry->current);
    entry->currentCoefficient(mRing, mLeadTerm.coeff);
    
    while (true) {
      ++entry->pos;
      if (entry->pos == entry->end) {
        mQueue.pop();
        entry->destroy(mRing);
        mPool.free(entry);
      } else {
        entry->computeCurrent(mRing);
        mQueue.decreaseTop(entry);
      }
      
      if (mQueue.empty())
        break;
      
      entry = mQueue.top();
      if (!mRing.monomialEQ(entry->current, mLeadTerm.monom))
        break;
      coefficient coeff;
      entry->currentCoefficient(mRing, coeff);
      mRing.coefficientAddTo(mLeadTerm.coeff, const_cast<const coefficient&>(coeff));
    }
  } while (mRing.coefficientIsZero(mLeadTerm.coeff));

  result = mLeadTerm;
  mLeadTermKnown = true;
  return true;
}

void TournamentReducer::removeLeadTerm()
{
  if (!mLeadTermKnown) {
    const_term dummy;
    leadTerm(dummy);
  }
  mLeadTermKnown = false;
}

void TournamentReducer::value(Poly &result)
{
  const_term t;
  while (leadTerm(t)) {
    result.appendTerm(t.coeff, t.monom);
    removeLeadTerm();
  }
  resetReducer();
}

void TournamentReducer::resetReducer()
{
  MonomialFree freeer(mRing);
  mQueue.forAll(freeer);
  mQueue.clear();
}

size_t TournamentReducer::getMemoryUse() const
{
  return
    TypicalReducer::getMemoryUse() +
    mQueue.getMemoryUse() +
    mPool.getMemoryUse();
}

void TournamentReducer::dump() const
{
}

MATHICGB_NAMESPACE_END
