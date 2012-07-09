// Copyright 2011 Michael E. Stillman

#ifndef _Reducer_h_
#define _Reducer_h_

#include "PolyRing.hpp"
#include "Poly.hpp"
#include <memtailor.h>
#include <memory>

class GroebnerBasis;
class PolyBasis;

/** Abstract base class for classes that allow reduction of polynomials.

todo: consider changing name of findLeadTerm to leadTerm.
*/
class Reducer {
public:
  enum ReducerType {
    Reducer_PolyHeap,
    Reducer_PolyGeoBucket,
    Reducer_Poly,
    Reducer_PolyHash,
    Reducer_BjarkeGeo, // uses hash table on front to remove duplicates
    Reducer_TournamentTree,
    Reducer_HashTourTree,

    Reducer_TourTree_NoDedup,
    Reducer_TourTree_Dedup,
    Reducer_TourTree_Hashed,
    Reducer_TourTree_NoDedup_Packed,
    Reducer_TourTree_Dedup_Packed,
    Reducer_TourTree_Hashed_Packed,

    Reducer_Heap_NoDedup,
    Reducer_Heap_Dedup,
    Reducer_Heap_Hashed,
    Reducer_Heap_NoDedup_Packed,
    Reducer_Heap_Dedup_Packed,
    Reducer_Heap_Hashed_Packed,

    Reducer_Geobucket_NoDedup,
    Reducer_Geobucket_Dedup,
    Reducer_Geobucket_Hashed,
    Reducer_Geobucket_NoDedup_Packed,
    Reducer_Geobucket_Dedup_Packed,
    Reducer_Geobucket_Hashed_Packed
  };

  struct Stats {
    Stats();

    // Number of signature reductions performed including singular reductions.
    unsigned long long reductions;

    // Number of regular reductions where the polynomial was not top
    // regular reducible.
    unsigned long long singularReductions;

    // Number of reductions to zero.
    unsigned long long zeroReductions;

    // Total number of steps across all reductions. Does not count detecting
    // singular reduction as a reduction step.
    unsigned long long steps;

    // Number of reduction steps in the reduction that had the most steps.
    unsigned long long maxSteps;
  };

  Reducer();
  virtual ~Reducer() {}

  static ReducerType reducerType(int typ);
  static void displayReducerTypes(std::ostream &o);

  Stats sigStats() const {return mSigStats;}
  Stats classicStats() const {return mClassicStats;}

  void reset();

  virtual std::string description() const = 0;
  virtual void insertTail(const_term multiplier, const Poly *f) = 0;
  virtual void insert(monomial multiplier, const Poly *f) = 0;

  virtual bool findLeadTerm(const_term &result) = 0;
  virtual void removeLeadTerm() = 0;

  virtual void dump() const {}

  virtual size_t getMemoryUse() const = 0;

  // Regular reduce multiple*basisElement in signature sig by the
  // basis elements in basis.
  //
  // Returns null (0) if multiple*basisElement is not regular top
  // reducible. This indicates a singular reduction.
  Poly* regularReduce(
    const_monomial sig,
    const_monomial multiple,
    size_t basisElement,
    const GroebnerBasis& basis);

  // Clasically reduces poly by the basis elements of basis. The reduction
  // is classic in that no signatures are taken into account.
  std::auto_ptr<Poly> classicReduce(const Poly& poly, const PolyBasis& basis);

  // Clasically reduces poly by the basis elements of basis, except that the
  // lead term is not reduced. The reduction is classic in that no signatures
  // are taken into account.
  std::auto_ptr<Poly> classicTailReduce(const Poly& poly, const PolyBasis& basis);

  // Clasically reduces the S-polynomial between a and b.
  std::auto_ptr<Poly> classicReduceSPoly
    (const Poly& a, const Poly& b, const PolyBasis& basis);

  static std::auto_ptr<Reducer> makeReducer
    (ReducerType t, PolyRing const& ring);
  static std::auto_ptr<Reducer> makeReducerNullOnUnknown
    (ReducerType t, PolyRing const& ring);

protected:
  std::auto_ptr<Poly> classicReduce(const PolyBasis& basis);
  std::auto_ptr<Poly> classicReduce
    (std::auto_ptr<Poly> partialResult, const PolyBasis& basis);

  virtual void resetReducer() = 0;

  size_t stats_maxsize;
  size_t stats_maxsize_live;
  unsigned long long stats_n_inserts;
  unsigned long long stats_n_compares;

  Stats mSigStats;
  Stats mClassicStats;
  memt::Arena mArena;
};

#if 0
template<typename Queue, bool Deduplicate> class ReducerPacked;  // *,(0,1),1, although Dup and DeDup have 
  // different Entry types...
template<typename Queue, bool Deduplicate> class ReducerNotPacked;  // *,(0,1),0
template<typename Queue> class ReducerHashedNotPacked;  // *,2,0
template<typename Queue> class ReducerHashedPacked;  // *,2,1
#endif

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
