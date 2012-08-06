#include "stdinc.h"
#include "PairTriangle.hpp"

#include <limits>
#include <stdexcept>
#include "FreeModuleOrder.hpp"

PairTriangle::PairTriangle(const FreeModuleOrder& order, const PolyRing& ring, size_t queueType):
  mColumnCount(0),
  mOrder(order),
  mRing(ring),
  mPairQueue(*this) {
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

namespace {
  template<class PairIterator>
  class IndexIterator {
  public:
	typedef typename PairIterator::iterator_category iterator_category;
	typedef typename PairIterator::value_type value_type;
	typedef typename PairIterator::difference_type difference_type;
	typedef typename PairIterator::pointer pointer;
	typedef typename PairIterator::reference reference;

	IndexIterator(PairIterator pairIterator): mIterator(pairIterator) {}
	IndexIterator& operator++() {++mIterator; return *this;}
	size_t operator*() const {return mIterator->i;}
	difference_type operator-(const IndexIterator<PairIterator>& it) const {
	  return mIterator - it.mIterator;
	}
	bool operator==(const IndexIterator<PairIterator>& it) const {
	  return mIterator == it.mIterator;
	}
	bool operator!=(const IndexIterator<PairIterator>& it) const {
	  return mIterator != it.mIterator;
	}

  private:
	PairIterator mIterator;
  };
}

void PairTriangle::endColumn() {
  mOrder.destructiveSort(mPrePairs);
  typedef IndexIterator<std::vector<PreSPair>::const_iterator> Iter;
  mPairQueue.addColumnDescending
	(Iter(mPrePairs.begin()), Iter(mPrePairs.end()));

  ++mColumnCount;
  ASSERT(mColumnCount == columnCount());
  for (std::vector<PreSPair>::iterator it = mPrePairs.begin();
	   it != mPrePairs.end(); ++it)
    mRing.freeMonomial(it->signature);
  mPrePairs.clear();
}

void PairTriangle::pop() {
  mPairQueue.pop();
}

size_t PairTriangle::getMemoryUse() const {
  return mPrePairs.capacity() * sizeof(mPrePairs.front()) +
	mPairQueue.getMemoryUse();
}

std::string PairTriangle::name() const {
  return "todo";
  // return mPairQueue.name();
}

std::pair<size_t, size_t> PairTriangle::topPair() const {
  return mPairQueue.topPair();
}

// Returns the minimal orderBy monomial along all pairs. This is the orderBy
// monomial of topPair().
const_monomial PairTriangle::topOrderBy() const {
  return mPairQueue.topPairData();
}
