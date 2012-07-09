#include "stdinc.h"
#include "PairTriangle.hpp"

#include <limits>
#include <stdexcept>
#include "FreeModuleOrder.hpp"

SPairGroup::SPairGroup(
  size_t fixedIndex,
  monomial signature,
  std::vector<PreSPair>& prePairs,
  memt::Arena& arena
):
  mSignature(signature),
  mFixedIndex(fixedIndex)
{
  const size_t prePairCount = prePairs.size();
  if (big()) {
    std::pair<BigIndex*, BigIndex*> range =
      arena.allocArrayNoCon<BigIndex>(prePairCount);
    for (size_t i = 0; i < prePairCount; ++i) {
      ASSERT(prePairs[i].i < mFixedIndex);
      ASSERT(prePairs[i].i <= std::numeric_limits<BigIndex>::max());
      range.first[i] = prePairs[i].i;
    }
    bigBegin = range.first;
    bigEnd = range.second;
  } else {
    std::pair<SmallIndex*, SmallIndex*> range =
      arena.allocArrayNoCon<SmallIndex>(prePairCount);
    for (size_t i = 0; i < prePairCount; ++i) {
      ASSERT(prePairs[i].i < mFixedIndex);
      ASSERT(prePairs[i].i <= std::numeric_limits<SmallIndex>::max());
      range.first[i] = prePairs[i].i;
    }
    smallBegin = range.first;
    smallEnd = range.second;
  }
  ASSERT(size() == prePairCount);
  ASSERT(empty() == (prePairCount == 0));
}

void SPairGroup::write(const PolyRing* R, std::ostream& out) {
  if (big()) {
    for (BigIndex* f = bigBegin; f != bigEnd; f++) {
      out << "  [" << fixedIndex() << " " << *f;
      if (f == bigBegin)
        R->monomialDisplay(out, signature());
      out << "]" << std::endl;
    }
  } else {
    for (SmallIndex* f = smallBegin; f != smallEnd; f++) {
      out << "  [" << fixedIndex() << " " << *f;
      if (f == smallBegin)
        R->monomialDisplay(out, signature());
      out << "]" << std::endl;
    }
  }
}

size_t SPairGroup::otherIndex() const {
  ASSERT(!empty());
  if (big())
    return *bigBegin;
  else
    return *smallBegin;
}

void SPairGroup::increment() {
  ASSERT(!empty());
  do {
    if (big())
      ++bigBegin;
    else
      ++smallBegin;
  } while (!empty() && otherIndex() == fixedIndex());
}

bool SPairGroup::empty() const {
  if (big())
    return bigBegin == bigEnd;
  else
    return smallBegin == smallEnd;
}

bool SPairGroup::big() const {
  return mFixedIndex >
    static_cast<size_t>(std::numeric_limits<SmallIndex>::max());
}

size_t SPairGroup::size() const {
  if (big())
    return bigEnd - bigBegin;
  else
    return smallEnd - smallBegin;
}

PairTriangle::PairTriangle(const FreeModuleOrder& order, const PolyRing& ring, size_t queueType):
  mUseSingletonGroups((queueType & 2) != 0),
  mColumnCount(0),
  mQueue(order.makeQueue((queueType & 1))),
  mOrder(order),
  mRing(ring) {
}

PairTriangle::~PairTriangle() {
}

size_t PairTriangle::size() const {
  size_t sum = 0;
  for (size_t i = 0; i < mGroups.size(); ++i)
    if (mGroups[i] != 0)
      sum += mGroups[i]->size();
  return sum;
}

void PairTriangle::beginColumn() {
  ASSERT(mPrePairs.empty());
  size_t const maxBigIndex = std::numeric_limits<BigIndex>::max();
  if (mColumnCount >= maxBigIndex)
    throw std::overflow_error
      ("Too large basis element index in constructing S-pairs.");
}

void PairTriangle::addPair(size_t index, monomial orderBy) {
  ASSERT(index < mColumnCount);
#ifdef DEBUG
  monomial tmp = mRing.allocMonomial();
  calculateOrderBy(mColumnCount, index, tmp);
  ASSERT(mRing.monomialEQ(tmp, orderBy));
  mRing.freeMonomial(tmp);
#endif

  PreSPair p;
  p.i = static_cast<BigIndex>(index);
  p.signature = orderBy;
  mPrePairs.push_back(p);
}

void PairTriangle::endColumn() {
  ++mColumnCount;
  if (mPrePairs.empty())
    return;

  size_t const newGen = mColumnCount - 1;
  if (mUseSingletonGroups) {
    size_t const size = mPrePairs.size();
    for (size_t i = 0; i < size; ++i) {
      mGroups.push_back(mArena.allocObjectNoCon<SPairGroup>());
      SPairGroup* p = mGroups.back();

      ASSERT(mPrePairTmp.empty());
      mPrePairTmp.push_back(mPrePairs[i]);
      new (p) SPairGroup
        (newGen, mPrePairTmp.front().signature, mPrePairTmp, mArena);
      mPrePairTmp.clear();

      calculateOrderBy(p->fixedIndex(), p->otherIndex(), p->signature());
      mQueue->push(p);
    }
  } else {
    mGroups.push_back(mArena.allocObjectNoCon<SPairGroup>());
    SPairGroup* p = mGroups.back();

    mOrder.destructiveSort(mPrePairs);

    // Take the very first one and insert it into the queue
    new (p) SPairGroup(newGen, mPrePairs.front().signature, mPrePairs, mArena);
    const size_t size = mPrePairs.size();
    for (size_t i = 1; i < size; ++i)
      mRing.freeMonomial(mPrePairs[i].signature);

    calculateOrderBy(p->fixedIndex(), p->otherIndex(), p->signature());
    mQueue->push(p);
  }
  mPrePairs.clear();
}

void PairTriangle::pop() {
  ASSERT(!empty());

  SPairGroup* topGroup = mQueue->top();
  ASSERT(topGroup != 0);

  while (true) {
    topGroup->increment();
    if (topGroup->empty()) {
      mQueue->pop();
      mRing.freeMonomial(topGroup->signature());
      topGroup->setSignature(0);
    } else {
      // Compute the signature of the next S-pair in topGroup.
      if (!calculateOrderBy
        (topGroup->fixedIndex(), topGroup->otherIndex(), topGroup->signature()))
        continue;
      mQueue->decreaseTop(topGroup);
    }
    break;
  }
}

size_t PairTriangle::getMemoryUse() const {
  return mArena.getMemoryUse() + mQueue->getMemoryUse();
}

std::string PairTriangle::name() const {
  return mQueue->getName() + (mUseSingletonGroups ? "" : "-triangle");
}

std::pair<size_t, size_t> PairTriangle::topPair() const {
  ASSERT(!mQueue->empty());
  SPairGroup* group = mQueue->top();
  ASSERT(group != 0);
  ASSERT(!group->empty());
  return std::make_pair(group->fixedIndex(), group->otherIndex());
}

// Returns the minimal orderBy monomial along all pairs. This is the orderBy
// monomial of topPair().
const_monomial PairTriangle::topOrderBy() const {
  ASSERT(!mQueue->empty());
  SPairGroup* group = mQueue->top();
  ASSERT(group != 0);
  ASSERT(!group->empty());
  return group->signature();
}
