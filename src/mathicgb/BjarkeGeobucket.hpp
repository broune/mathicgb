// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_BJARKE_GEOBUCKET_GUARD
#define MATHICGB_BJARKE_GEOBUCKET_GUARD

#include "Reducer.hpp"
#include "PolyHashTable.hpp"
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

template<
  bool TrackFront,
  bool MinBucketBinarySearch,
  bool Deduplicate,
  bool Premerge,
  bool CollectMax,
  int BucketStorage,
  size_t InsertFactor = 1>
class GeoConfiguration {
public:
  GeoConfiguration(
    const PolyRing *R0,
    size_t geoBase,
    size_t minBucketSize
  ):
    R(R0), geoBase(geoBase), minBucketSize(minBucketSize), _comparisons(0) {}

  const PolyRing *R;
  size_t geoBase;
  size_t minBucketSize;

  typedef PolyHashTable::node * Entry;

  typedef bool CompareResult;

  CompareResult compare(const Entry& a, const Entry& b) const {
    ++_comparisons;
    return R->monomialLT(a->monom, b->monom);
  }
  bool cmpLessThan(CompareResult r) const {return r;}

  static const bool supportDeduplication = Deduplicate;
  bool cmpEqual(CompareResult r) const {return r;} // NOT USED IN OUR CASE HERRE!
  Entry deduplicate(const Entry& a, const Entry& /* b */) const {
    MATHICGB_ASSERT(false);
    return a;
  }

  size_t getComparisons() const {return _comparisons;}
  void resetComparisons() const {_comparisons = 0;}

  static const bool minBucketBinarySearch = MinBucketBinarySearch;
  static const bool trackFront = TrackFront;
  static const bool premerge = Premerge;
  static const bool collectMax = CollectMax;
  static const mic::GeobucketBucketStorage bucketStorage =
    (mic::GeobucketBucketStorage)BucketStorage;
  static const size_t insertFactor = InsertFactor;

private:
  mutable size_t _comparisons;
};

class BjarkeGeobucket : public Reducer {
public:
  BjarkeGeobucket(const PolyRing *R);
  ~BjarkeGeobucket();

  virtual std::string description() const { return "bjarke geo buckets"; }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  bool findLeadTerm(const_term &result);
  void removeLeadTerm();

  void value(Poly &result); // keep extracting lead term until done

  size_t getMemoryUse() const;

  void dump() const; // Used for debugging

protected:
  void resetReducer();

private:
  void insert(const_term multiplier, Poly::iterator first, Poly::iterator last);

  typedef PolyHashTable::MonomialArray HashPoly;

  const PolyRing *R_;
  PolyHashTable H_;

  size_t mNodeCount;  // number of (distinct) monomials in G_

  typedef GeoConfiguration<true,true,false,false,false,1,1> Configuration;
  mic::Geobucket< Configuration > G_;
};

MATHICGB_NAMESPACE_END
#endif
