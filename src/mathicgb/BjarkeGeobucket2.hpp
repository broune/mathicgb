// Copyright 2011 Michael E. Stillman

#ifndef _bjarkeGeoBucket_h_
#define _bjarkeGeoBucket_h_

#include <mathic.h>
#include "TypicalReducer.hpp"
#include "PolyHashTable.hpp"

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

  static const bool supportDeduplication = false;
  bool cmpEqual(CompareResult r) const {ASSERT(false);return r;} // NOT USED IN OUR CASE HERRE!
  Entry deduplicate(const Entry& a, const Entry& /* b */) const {ASSERT(false); return a;}

  size_t getComparisons() const {return _comparisons;}
  void resetComparisons() const {_comparisons = 0;}

  static const bool minBucketBinarySearch = true; // MinBucketBinarySearch;
  static const bool trackFront = true; //TrackFront;
  static const bool premerge = false;
  static const bool collectMax = false;
  static const mic::GeobucketBucketStorage bucketStorage = static_cast<mic::GeobucketBucketStorage>(1);
  static const size_t insertFactor = 1;

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

class BjarkeGeobucket2 : public TypicalReducer {
public:
  BjarkeGeobucket2(const PolyRing *R);

  virtual std::string description() const { return "bjarke geo buckets"; }

  void insertTail(const_term multiplier, const Poly *f);
  void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term &result);
  virtual void removeLeadTerm();

  void value(Poly &result); // keep extracting lead term until done

  virtual size_t getMemoryUse() const;

  void dump() const; // Used for debugging

protected:
  void resetReducer();

private:
  typedef mic::HashTable<BjarkeGeobucket2Configuration>::Handle node;
  typedef PolyHashTable::MonomialArray HashPoly;

  void insert(Poly::const_iterator first, 
              Poly::const_iterator last,
              std::vector<node*> &result);

  void insert(const_term multiplier, Poly::iterator first, Poly::iterator last);


  const PolyRing &mRing;
  PolyHashTable mHashTableOLD;
  mic::HashTable<BjarkeGeobucket2Configuration> mHashTable;
  mic::Geobucket< GeoConfiguration > mHeap;
};

#endif

// Local Variables:
// indent-tabs-mode: nil
// compile-command: "make -C $MATHIC/mathicgb "
// End:
