#include "stdinc.h"
#include "SPairs.hpp"

#include "GroebnerBasis.hpp"

#include <iostream>

SPairs::SPairs(const PolyBasis& basis, size_t queueType):
  mTri(basis, queueType),
  mBasis(basis),
  mRing(basis.ring()) {}

std::pair<size_t, size_t> SPairs::pop() {
  // Must call addPairs for new elements before popping.
  ASSERT(mEliminated.columnCount() == mBasis.size());

  while (!mTri.empty()) {
    std::pair<size_t, size_t> p;
    p = mTri.topPair();
    if (mBasis.retired(p.first) || mBasis.retired(p.second)) {
      mTri.pop();
      continue;
    }
    const_monomial lcm = mTri.topOrderBy();
    ASSERT(mRing.monomialIsLeastCommonMultiple
      (mBasis.leadMonomial(p.first),
      mBasis.leadMonomial(p.second), lcm));
    // Can't pop before done with lcm as popping overwrites lcm.
    if (!advancedBuchbergerLcmCriterion(p.first, p.second, lcm)) {
      mTri.pop();
      mEliminated.setBit(p.first, p.second, true);
      return p;
    }
    mTri.pop();
  }
  return std::make_pair(static_cast<size_t>(-1), static_cast<size_t>(-1));
}

size_t SPairs::size() const {
  return mTri.size();
}

namespace {
  // Records multiples of a basis element.
  // Used in addPairs().
  class RecordIndexes : public DivisorLookup::EntryOutput {
  public:
    RecordIndexes(
      size_t newGen,
      BitTriangle& eliminated,
      std::vector<size_t>& indexes
    ):
      mNewGen(newGen),
      mEliminated(eliminated),
      mIndexes(indexes) {
      mIndexes.clear();
    }
    virtual bool proceed(size_t index) {
      if (index == mNewGen)
        return true;
      mIndexes.push_back(index);

      // The S-pair (newGen, *it) corresponds to reducing the non-minimal
      // basis element *it. The S-polynomial corresponds to the first
      // step of that reduction. We tell the caller to reduce *it, so we
      // get to assume that this S-pair can be eliminated. This is important
      // because it sometimes allows us to eliminate an S-pair (newGen, x)
      // when (*it, x) has already been eliminated. We want to make use of
      // this opportunity before removing all the information about *it.
      mEliminated.setBit(mNewGen, index, true);
      return true;
    }
  private:
    size_t const mNewGen;
    BitTriangle& mEliminated;
    std::vector<size_t>& mIndexes;
  };
}

void SPairs::addPairsAssumeAutoReduce(
  size_t newGen,
  std::vector<size_t>& toRetireAndReduce
) {
  ASSERT(mTri.columnCount() == newGen);

  ASSERT(newGen < mBasis.size());
  ASSERT(!mBasis.retired(newGen));

  while (mEliminated.columnCount() < mBasis.size()) {
    if (mUseBuchbergerLcmHitCache) {
      ASSERT(mEliminated.columnCount() == mBuchbergerLcmHitCache.size());
      mBuchbergerLcmHitCache.push_back(0);
    }
    mEliminated.addColumn();
  }

  RecordIndexes indexes(newGen, mEliminated, toRetireAndReduce);
  mBasis.divisorLookup().multiples(mBasis.leadMonomial(newGen), indexes);
  addPairs(newGen);
}

void SPairs::addPairs(size_t newGen) {
  // Must call addPairs with newGen parameter in the sequence 0, 1, ...
  // newGen could be implicitly picked up from mTri.columnCount(), but
  // doing it this way ensures that what happens is what the client thinks
  // is happening and offers an ASSERT to inform mistaken client code.
  ASSERT(mTri.columnCount() == newGen);

  ASSERT(newGen < mBasis.size());
  ASSERT(!mBasis.retired(newGen));

  while (mEliminated.columnCount() < mBasis.size()) {
    if (mUseBuchbergerLcmHitCache) {
      ASSERT(mEliminated.columnCount() == mBuchbergerLcmHitCache.size());
      mBuchbergerLcmHitCache.push_back(0);
    }
    mEliminated.addColumn();
  }

  mTri.beginColumn();
  const_monomial const newLead = mBasis.leadMonomial(newGen);
  monomial lcm = mBasis.ring().allocMonomial();
  for (size_t oldGen = 0; oldGen < newGen; ++oldGen) {
    if (mBasis.retired(oldGen))
      continue;
    const_monomial const oldLead = mBasis.leadMonomial(oldGen);
    if (mRing.monomialRelativelyPrime(newLead, oldLead)) {
      ++mStats.relativelyPrimeHits;
      mEliminated.setBit(newGen, oldGen, true);
      continue;
    }
    mRing.monomialLeastCommonMultipleNoWeights(newLead, oldLead, lcm);
    if (simpleBuchbergerLcmCriterion(newGen, oldGen, lcm)) {
      mEliminated.setBit(newGen, oldGen, true);
      continue;
    }
    mRing.setWeightsOnly(lcm);
    mTri.addPair(oldGen, lcm);
    lcm = mBasis.ring().allocMonomial();
  }
  mBasis.ring().freeMonomial(lcm);
  mTri.endColumn();
}

size_t SPairs::getMemoryUse() const {
  return mTri.getMemoryUse();
}

SPairs::ClassicPairTriangle::ClassicPairTriangle(const PolyBasis& basis, size_t queueType):
  PairTriangle(basis.order(), basis.ring(), queueType),
  mBasis(basis) {}

bool SPairs::simpleBuchbergerLcmCriterion
  (size_t a, size_t b, const_monomial lcmAB) const
{
  ASSERT(a < mBasis.size());
  ASSERT(b < mBasis.size());
  ASSERT(a != b);
  ASSERT(!mBasis.retired(a));
  ASSERT(!mBasis.retired(b));
  ASSERT(mRing.monomialIsLeastCommonMultipleNoWeights
    (mBasis.leadMonomial(a), mBasis.leadMonomial(b), lcmAB));
  ASSERT(mEliminated.columnCount() == mBasis.size());

  class Criterion : public DivisorLookup::EntryOutput {
  public:
    Criterion(size_t a, size_t b, const_monomial lcmAB, const SPairs& sPairs):
      mA(a), mB(b),
      mLcmAB(lcmAB),
      mSPairs(sPairs),
      mRing(sPairs.ring()),
      mBasis(sPairs.basis()),
      mHit(static_cast<size_t>(-1)),
      mAlmostApplies(false) {}

    virtual bool proceed(size_t index) {
      ASSERT(index < mBasis.size());
      ASSERT(!applies()); // should have stopped search in this case
      ASSERT(mRing.monomialIsDivisibleBy(mLcmAB, mBasis.leadMonomial(index)));
      if (!mBasis.leadMinimal(index))
        return true;
      if (index == mA || index == mB)
        return true;
      mAlmostApplies = true;

      // check lcm(a,index) != lcm(a,b) <=>
      // exists i such that max(a[i], c[i]) != max(a[i],b[i]) <=>
      // exists i such that b[i] > a[i] && b[i] > c[i]
      const_monomial leadA = mBasis.leadMonomial(mA);
      const_monomial leadB = mBasis.leadMonomial(mB);
      const_monomial leadC = mBasis.leadMonomial(index);
      if (!mSPairs.eliminated(index, mA) &&
        !mRing.monomialHasStrictlyLargerExponent(leadB, leadC, leadA))
        return true;

      // check lcm(b,index) != lcm(a,b)
      if (!mSPairs.eliminated(index, mB) &&
        !mRing.monomialHasStrictlyLargerExponent(leadA, leadC, leadB))
        return true;

      mHit = index;
      return false; // stop search
    }

    const_monomial lcmAB() const {return mLcmAB;}
    bool almostApplies() const {return mAlmostApplies;}
    bool applies() const {return mHit != static_cast<size_t>(-1);}
    size_t hit() const {return mHit;}

  private:
    size_t mA;
    size_t mB;
    const_monomial mLcmAB;
    const SPairs& mSPairs;
    const PolyRing& mRing;
    const PolyBasis& mBasis;
    size_t mHit; // the divisor that made the criterion apply
    bool mAlmostApplies; // applies ignoring lcm(a,b)=lcm(a,c) complication
  };

  bool applies = false;
  bool almostApplies = false;
  {
    Criterion criterion(a, b, lcmAB, *this);
    if (mUseBuchbergerLcmHitCache) {
      // Check cacheB first since when I tried this there was a higher hit rate
      // for cacheB than cacheA. Might not be a persistent phenomenon, but
      // there's no downside to trying out cacheB first so I'm going for that.
      //
      // I update the cache if the second check is a hit but not if the first
      // check is a hit. In the one test I did, the worst hit rate was from
      // updating the cache every time, the second best hit rate was from
      // not updating the cache (from cache hits) and the best hit rate was
      // from doing this.
      //
      // The idea is that when the first cache check is a hit,
      // the second cache member might have been a hit too, and updating it
      // might replace a high hit rate element with a low hit rate element,
      // which would be bad. When the second cache check is a hit, we know
      // that the first one wasn't (or we would have taken an early exit),
      // so we have reason to suspect that the first cache element is not
      // a high hit rate element. So it should be better to replace it.
      // That idea seems to be right since it worked better in the one
      // test I did.
      size_t cacheB = mBuchbergerLcmHitCache[b];
      if (!applies && !mBasis.retired(cacheB) &&
        mRing.monomialIsDivisibleBy
          (criterion.lcmAB(), mBasis.leadMonomial(cacheB)))
        applies = !criterion.Criterion::proceed(cacheB);

      size_t cacheA = mBuchbergerLcmHitCache[a];
      if (!applies && !mBasis.retired(cacheA) &&
        mRing.monomialIsDivisibleBy
          (criterion.lcmAB(), mBasis.leadMonomial(cacheA))) {
        applies = !criterion.Criterion::proceed(cacheA);
        if (applies)
          mBuchbergerLcmHitCache[b] = cacheA;
      }
    }
    if (applies)
      {
	if (mStats.late)
	  ++mStats.buchbergerLcmCacheHitsLate;
	else
	  ++mStats.buchbergerLcmCacheHits;
      }
    else {
      ASSERT(!criterion.applies());
      mBasis.divisorLookup().divisors(criterion.lcmAB(), criterion);
      applies = criterion.applies();

      if (mUseBuchbergerLcmHitCache && applies) {
        ASSERT(criterion.hit() < mBasis.size());
        mBuchbergerLcmHitCache[a] = criterion.hit();
        mBuchbergerLcmHitCache[b] = criterion.hit();
      }
    }
    if (!applies)
      almostApplies = criterion.almostApplies();
  }

  if (applies)
    {
      if (mStats.late)
	++mStats.buchbergerLcmSimpleHitsLate;
      else
	++mStats.buchbergerLcmSimpleHits;
    }

  ASSERT(applies == simpleBuchbergerLcmCriterionSlow(a, b));
  return applies;
}

bool SPairs::simpleBuchbergerLcmCriterionSlow(size_t a, size_t b) const {
  ASSERT(a < mBasis.size());
  ASSERT(b < mBasis.size());
  ASSERT(a != b);
  ASSERT(!mBasis.retired(a));
  ASSERT(!mBasis.retired(b));
  ASSERT(mEliminated.columnCount() == mBasis.size());

  // todo: use iterators
  monomial lcmAB = mRing.allocMonomial();
  monomial lcm = mRing.allocMonomial();
  mRing.monomialLeastCommonMultiple
    (mBasis.leadMonomial(a), mBasis.leadMonomial(b), lcmAB);
  size_t stop = mBasis.size();
  size_t i = 0;
  for (; i < stop; ++i) {
    if (mBasis.retired(i) || !mBasis.leadMinimal(i))
      continue;
    if (!mRing.monomialIsDivisibleBy(lcmAB, mBasis.leadMonomial(i)))
      continue;
    if (i == a || i == b)
      continue;
    if (!eliminated(i, a)) {
      mRing.monomialLeastCommonMultiple
        (mBasis.leadMonomial(a), mBasis.leadMonomial(i), lcm);
      if (mRing.monomialEQ(lcmAB, lcm))
        continue;
    }

    if (!eliminated(i, b)) {
      mRing.monomialLeastCommonMultiple
        (mBasis.leadMonomial(b), mBasis.leadMonomial(i), lcm);
      if (mRing.monomialEQ(lcmAB, lcm))
        continue;
    }
    break;
  }
  mRing.freeMonomial(lcmAB);
  mRing.freeMonomial(lcm);
  return i != stop;
}

bool SPairs::advancedBuchbergerLcmCriterion
  (size_t a, size_t b, const_monomial lcmAB) const
{
  ASSERT(a != b);
  ASSERT(mEliminated.columnCount() == mBasis.size());
  ASSERT(mRing.monomialIsLeastCommonMultipleNoWeights
    (mBasis.leadMonomial(a), mBasis.leadMonomial(b), lcmAB));

  mStats.late = true;
  if (simpleBuchbergerLcmCriterion(a, b, lcmAB)) {
    mStats.late = false;
    ASSERT(advancedBuchbergerLcmCriterionSlow(a, b));
    return true;
  }
  mStats.late = false;

  // *** Build the graph vertices
  // graph contains pairs (index, state). index is the index of a basis
  // that is a node in G. state indicates which of a and b that the node
  // in question is so far known to be connected to.

  typedef std::vector<std::pair<size_t, Connection> > Graph;
  class GraphBuilder : public DivisorLookup::EntryOutput {
  public:
    GraphBuilder(Graph& graph): mGraph(graph) {graph.clear();}
    virtual bool proceed(size_t index) {
      mGraph.push_back(std::make_pair(index, NotConnected));
      return true;
    }
  private:
    Graph& mGraph;
  };
  Graph& graph = mAdvancedBuchbergerLcmCriterionGraph;
  graph.clear();
  GraphBuilder builder(graph);
  mBasis.divisorLookup().divisors(lcmAB, builder);

  if (graph.size() <= 3) {
    // For the graph approach to be better than the simpler approach of
    // considering triples, there has to be more than 3 nodes in the graph.
    ASSERT(!advancedBuchbergerLcmCriterionSlow(a, b));
    return false;
  }

  // *** Set up todo with a and b
  // todo points to elements (nodes) of graph to process.
  std::vector<Graph::iterator> todo;
  Graph::iterator const graphEnd = graph.end();
  for (Graph::iterator it = graph.begin(); it != graphEnd; ++it) {
    if (it->first == a)
      it->second = ConnectedA;
    else if (it->first == b)
      it->second = ConnectedB;
    else
      continue;
    todo.push_back(it);
  }

  // *** Follow edges in the graph
  // We stop as soon as we find a node that is connected to both a and b,
  // since then a and b are connected so that the criterion applies.
  bool applies = false;
  while (!applies && !todo.empty()) {
    size_t const currentIndex = todo.back()->first;
    Connection const currentConnect = todo.back()->second;
    ASSERT(currentConnect != NotConnected);
    todo.pop_back();

    // loop through all potential edges (currentIndex, otherIndex)
    const_monomial const currentLead = mBasis.leadMonomial(currentIndex);
    for (Graph::iterator other = graph.begin(); other != graphEnd; ++other) {
      Connection const otherConnect = other->second;
      if (currentConnect == otherConnect)
        continue;
      size_t const otherIndex = other->first;
      ASSERT(otherIndex != currentIndex);

      const_monomial const otherLead = mBasis.leadMonomial(otherIndex);
      // Note that
      //  lcm(c,d) != lcmAB <=>
      //  exists i such that max(c[i], d[i]) < lcmAB[i] <=>
      //  exists i such that lcmAB[i] > c[i] && lcmAB[i] > d[i]
      if (!eliminated(currentIndex, otherIndex) &&
        !mRing.monomialHasStrictlyLargerExponent(lcmAB, currentLead, otherLead))
      {
        continue; // not an edge in G
      }

      if (otherConnect == NotConnected) {
        other->second = currentConnect;
        todo.push_back(other);
      } else {
        // At this point we have found an edge between a node connected to
        // a and a node connected to b. So a and b are connected.
        ASSERT(currentConnect != otherConnect);
        applies = true;
        break;
      }
    }
  }

  if (applies)
    ++mStats.buchbergerLcmAdvancedHits;

  //  if (graph.size() >= 10)
  //    std::cout << "[adv size=" << graph.size() << " result= " << applies << std::endl;

  ASSERT(applies == advancedBuchbergerLcmCriterionSlow(a, b));
  return applies;
}

bool SPairs::advancedBuchbergerLcmCriterionSlow(size_t a, size_t b) const {
  ASSERT(a != b);
  ASSERT(mEliminated.columnCount() == mBasis.size());

  monomial lcmAB = mRing.allocMonomial();
  monomial lcm = mRing.allocMonomial();
  mRing.monomialLeastCommonMultiple
    (mBasis.leadMonomial(a), mBasis.leadMonomial(b), lcmAB);
  size_t stop = mBasis.size();

  // *** Build the graph vertices
  // graph contains pairs (index, state). index is the index of a basis
  // that is a node in G. state indicates which of a and b that the node
  // in question is so far known to be connected to.
  std::vector<std::pair<size_t, Connection> > graph;
  std::vector<size_t> todo; // indexes into graph to process.
  for (size_t i = 0; i != stop; ++i) {
    if (mBasis.retired(i))
      continue;
    if (!mRing.monomialIsDivisibleBy(lcmAB, mBasis.leadMonomial(i)))
      continue;
    Connection con = NotConnected;
    if (i == a) {
      con = ConnectedA;
      todo.push_back(graph.size());
    } else if (i == b) {
      con = ConnectedB;
      todo.push_back(graph.size());
    }
    graph.push_back(std::make_pair(i, con));
  }

  // *** Follow edges in the graph
  // We stop as soon as we find a node that is connected to both a and b,
  // since then a and b are connected so that the criterion applies.
  bool applies = false;
  while (!applies && !todo.empty()) {
    ASSERT(todo.size() <= graph.size());
    std::pair<size_t, Connection> const node = graph[todo.back()];
    todo.pop_back();
    ASSERT(node.second != NotConnected);

    // loop through all potential edges (node.first, i)
    const_monomial leadNode = mBasis.leadMonomial(node.first);
    for (size_t i = 0; i < graph.size(); ++i) {
      if (node.second == graph[i].second)
        continue;
      ASSERT(graph[i].first != node.first);
      size_t const other = graph[i].first;

      const_monomial const leadOther = mBasis.leadMonomial(other);
      mRing.monomialLeastCommonMultiple(leadNode, leadOther, lcm);
      if (!eliminated(node.first, other) && mRing.monomialEQ(lcm, lcmAB))
        continue; // not an edge in G
      
      if (graph[i].second == NotConnected) {
        graph[i].second = node.second;
        todo.push_back(i);
      } else {
        // At this point we have found an edge between something a node to
        // a and a node connected to b. So a and b are connected.
        ASSERT(graph[i].second != node.second)
        applies = true;
        break;
      }
    }
  }
  mRing.freeMonomial(lcmAB);
  mRing.freeMonomial(lcm);

  ASSERT(applies || !simpleBuchbergerLcmCriterionSlow(a, b));
  return applies;
}

SPairs::Stats SPairs::stats() const {
  size_t const columnCount = mTri.columnCount();
  mStats.sPairsConsidered = columnCount * (columnCount - 1) / 2;
  return mStats;
}

std::string SPairs::name() const {
  return mTri.name();
}

bool SPairs::ClassicPairTriangle::calculateOrderBy(
  size_t a,
  size_t b,
  monomial orderBy
) const {
  ASSERT(!orderBy.isNull());
  ASSERT(a != b);
  ASSERT(a < mBasis.size());
  ASSERT(b < mBasis.size());
  if (mBasis.retired(a) || mBasis.retired(b))
    return false;
  const_monomial const leadA = mBasis.leadMonomial(a);
  const_monomial const leadB = mBasis.leadMonomial(b);
  mBasis.ring().monomialLeastCommonMultiple(leadA, leadB, orderBy);
  return true;
}
