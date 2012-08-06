#ifndef _s_pairs_h_
#define _s_pairs_h_

#include "PairTriangle.hpp"
#include <utility>
#include <mathic.h>

class PolyBasis;

// Stores the set of pending S-pairs for use in the classic Buchberger
// algorithm. Also eliminates useless S-pairs and orders the S-pairs.
class SPairs {
public:
  SPairs(const PolyBasis& basis, size_t queueType);

  // Returns the number of S-pairs in the data structure.
  size_t pairCount() const {mTri.pairCount();}

  // Returns true if there all S-pairs have been eliminated or popped.
  bool empty() const {return mTri.empty();}

  // Removes the minimal S-pair from the data structure and return it.
  // The S-polynomial of that pair is assumed to reduce to zero, either
  // because it already does, or because it did not reduce to zero and then
  // that caused the addition of another basis element. If this assumption
  // is broken, S-pair elimination may give incorrect results.
  //
  // Returns the pair (invalid,invalid) if there are no S-pairs left to
  // return, where invalid is static_cast<size_t>(-1). This can happen even
  // if empty() returned false prior to calling pop(), since the S-pairs in
  // the queue may have been found to be useless.
  std::pair<size_t, size_t> pop();

  // Add the pairs (index,a) to the data structure for those a such that
  // a < index. Some of those pairs may be eliminated if they can be proven
  // to be useless. index must be a valid index of a basis element
  // and must only be called once per basis element and must be 1 more
  // than the value of index for the previous call to addPairs, starting
  // at zero for the first call.
  void addPairs(size_t index);

  // As addPairs, but assuming auto-reduction of the basis will happen.
  // This method assumes that if lead(index) divides lead(x) for a basis
  // element x, then x will be retired from the basis and reduced. toReduce
  // will contain those indices x.
  void addPairsAssumeAutoReduce(size_t index, std::vector<size_t>& toRetireAndReduce);

  // Returns true if the S-pair (a,b) is known to be useless. Even if the
  // S-pair is not useless now, it will become so later. At the latest, an
  // S-pair becomes useless when its S-polynomial has been reduced to zero.
  bool eliminated(size_t a, size_t b) const {
    return mEliminated.bitUnordered(a, b);
  }

  const PolyRing& ring() const {return mRing;}
  const PolyBasis& basis() const {return mBasis;}

  size_t getMemoryUse() const;

  struct Stats {
    Stats():
      sPairsConsidered(0),
      relativelyPrimeHits(0),
      buchbergerLcmSimpleHits(0),
      buchbergerLcmAdvancedHits(0),
      buchbergerLcmCacheHits(0),
      late(false),
      buchbergerLcmSimpleHitsLate(0),
      buchbergerLcmCacheHitsLate(0)
    {}

    unsigned long long sPairsConsidered;
    unsigned long long relativelyPrimeHits;
    unsigned long long buchbergerLcmSimpleHits;
    unsigned long long buchbergerLcmAdvancedHits;
    unsigned long long buchbergerLcmCacheHits;
    bool late;  // if set to true then simpleBuchbergerLcmCriterion sets the following 2 instead:
    unsigned long long buchbergerLcmSimpleHitsLate;
    unsigned long long buchbergerLcmCacheHitsLate;
  };
  Stats stats() const;

  std::string name() const;

private:
  // Returns true if Buchberger's second criterion for eliminating useless
  // S-pairs applies to the pair (a,b). Let
  //   l(a,b) = lcm(lead(a), lead(b)).
  // The criterion says that if there is some other basis element c such that
  //   lead(c)|l(a,b)
  // and
  //   l(a,c) reduces to zero, and
  //   l(b,c) reduces to zero
  // then (a,b) will reduce to zero (using classic non-signature reduction).
  //
  // This criterion is less straight forward to apply in case for example
  //   l(a,b) = l(a,c) = l(b,c)
  // since then there is the potential to erroneously eliminate all the three
  // pairs among a,b,c on the assumption that the other two pairs will reduce
  // to zero. In such cases, we eliminate the pair with the lowest indexes.
  // This allows removing generators that get non-minimal lead term without
  // problems.
  bool simpleBuchbergerLcmCriterion
    (size_t a, size_t b, const_monomial lcmAB) const;

  // As the non-slow version, but uses simpler and slower code.
  bool simpleBuchbergerLcmCriterionSlow(size_t a, size_t b) const;

  // Improves on Buchberger's second criterion by using connection in a graph
  // to determine if an S-pair can be eliminated. This can eliminate some pairs
  // that cannot be eliminated by looking at any one triple of generators.
  //
  // The algorithm is based on considering an undirected graph G.
  // Each vertex of G represents a basis element whose lead monomial divides
  // lcmAB. There is an edge (c,d) if lcm(c,d) != lcm(a,b) or if (c,d) has
  // been eliminated. It is a theorem that if there is a path from a to b
  // in G then (a,b) is a useless S-pair and so can be eliminated.
  bool advancedBuchbergerLcmCriterion
    (size_t a, size_t b, const_monomial lcmAB) const;

  // As the non-slow version, but uses simpler and slower code.
  bool advancedBuchbergerLcmCriterionSlow(size_t a, size_t b) const;

  class ClassicPairTriangle : public PairTriangle {
  public:
    ClassicPairTriangle(const PolyBasis& basis, size_t queueType);
  protected:
    virtual bool calculateOrderBy(size_t a, size_t b, monomial orderBy) const;
  private:
    const PolyBasis& mBasis;
  };
  ClassicPairTriangle mTri;

  // The bit at (i,j) is set to true if it is known that the S-pair between
  // basis element i and j does not have to be reduced. This can be due to a
  // useless S-pair criterion eliminating that pair, or it can be because the
  // S-polynomial of that pair has already been reduced.
  mathic::BitTriangle mEliminated;
  const PolyBasis& mBasis;
  const PolyRing& mRing;
  mutable Stats mStats;

  static const bool mUseBuchbergerLcmHitCache = true;
  mutable std::vector<size_t> mBuchbergerLcmHitCache;

  enum Connection { // used in advancedBuchbergerLcmCriterion().
    NotConnected, // not known to be connected to a or b
    ConnectedA, // connected to a
    ConnectedB // connected to b
  };
  // Variable used only inside advancedBuchbergerLcmCriterion().
  mutable std::vector<std::pair<size_t, Connection> >
    mAdvancedBuchbergerLcmCriterionGraph;
};

#endif
