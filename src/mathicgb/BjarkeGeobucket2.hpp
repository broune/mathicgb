// Copyright 2011 Michael E. Stillman

#ifndef _bjarkeGeoBucket_h_
#define _bjarkeGeoBucket_h_

#include <mathic.h>
#include "Reducer.hpp"
#include "PolyHashTable.hpp"

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
                   const PolyRing &ring,
                   size_t geoBase,
                   size_t minBucketSize):
    mRing(ring), geoBase(geoBase), minBucketSize(minBucketSize), _comparisons(0) {}

  const PolyRing &mRing;
  size_t geoBase;
  size_t minBucketSize;

  typedef PolyHashTable::node * Entry;

  typedef bool CompareResult;

  CompareResult compare(const Entry& a, const Entry& b) const {
    ++_comparisons;
    return mRing.monomialLT(a->monom, b->monom);
  }
  bool cmpLessThan(CompareResult r) const {return r;}

  static const bool supportDeduplication = Deduplicate;
  bool cmpEqual(CompareResult r) const {return r;} // NOT USED IN OUR CASE HERRE!
  Entry deduplicate(const Entry& a, const Entry& /* b */) const {ASSERT(false); return a;}

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

class BjarkeGeobucket2Configuration
{
public:
  typedef const_monomial Key;
  typedef coefficient Value;

  BjarkeGeobucket2Configuration(const PolyRing &ring) : mRing(ring) {}

  size_t hash(Key k) {return mRing.monomialHashValue(k);}

  bool keysEqual(Key k1, Key k2) {
    return ((mRing.monomialHashValue(k1) == mRing.monomialHashValue(k2)) 
            && mRing.monomialEQ(k1, k2));
  }

  void combine(Value &a, const Value &b) {mRing.coefficientAddTo(a,b);}
private:
  const PolyRing &mRing;
};

class BjarkeGeobucket2 : public Reducer {
public:
  BjarkeGeobucket2(const PolyRing *R);
  ~BjarkeGeobucket2();

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
  typedef GeoConfiguration<true,true,false,false,false,1,1> Configuration;

  size_t mNodeCount;  // number of (distinct) monomials in mHeap

  const PolyRing &mRing;
  PolyHashTable mHashTableOLD;
  mic::Geobucket< Configuration > mHeap;
};

#endif

// Local Variables:
// indent-tabs-mode: nil
// compile-command: "make -C $MATHIC/mathicgb "
// End:
