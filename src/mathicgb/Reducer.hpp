// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_REDUCER_GUARD
#define MATHICGB_REDUCER_GUARD

#include "PolyRing.hpp"
#include "Poly.hpp"
#include <memtailor.h>
#include <memory>

MATHICGB_NAMESPACE_BEGIN

class SigPolyBasis;
class PolyBasis;

/** Abstract base class for classes that allow reduction of polynomials.

todo: consider changing name of findLeadTerm to leadTerm.
*/
class Reducer {
public:
  virtual ~Reducer();

  /// Returns the preferred number of reductions to do at a time. A classic
  /// serial reducer will have no particular benefit from doing more than
  /// one reduction at a time, so it should say that. A matrix-based or
  /// parallel reducer will have benefit from being presented with
  /// larger sets of reductions at a time.
  virtual size_t preferredSetSize() const = 0;

  // ***** Methods that do reduction

  /** Clasically reduces poly by the basis elements of basis. The reduction
    is classic in that no signatures are taken into account. */
  virtual std::unique_ptr<Poly> classicReduce
  (const Poly& poly, const PolyBasis& basis) = 0;

  /** Clasically reduces poly by the basis elements of basis, except that the
   lead term is not reduced. The reduction is classic in that no signatures
   are taken into account. */
  virtual std::unique_ptr<Poly> classicTailReduce
  (const Poly& poly, const PolyBasis& basis) = 0;

  /** Clasically reduces the S-polynomial between a and b. */
  virtual std::unique_ptr<Poly> classicReduceSPoly
  (const Poly& a, const Poly& b, const PolyBasis& basis) = 0;

  /** Clasically reduces the S-polynomial of these pairs. May or may
      not also interreduce these to some extent. Polynomials that are
      reduced to zero are not put into reducedOut. */
  virtual void classicReduceSPolySet
  (std::vector<std::pair<size_t, size_t> >& spairs,
   const PolyBasis& basis,
   std::vector<std::unique_ptr<Poly> >& reducedOut) = 0;

  /** Clasically reduces the passed-in polynomials of these pairs. May
      or may not also interreduce these to some extent. Polynomials
      that are reduced to zero are not put into reducedOut. */
  virtual void classicReducePolySet
  (const std::vector<std::unique_ptr<Poly> >& polys,
   const PolyBasis& basis,
   std::vector<std::unique_ptr<Poly> >& reducedOut) = 0;

  /** Regular reduce multiple*basisElement in signature sig by the
    basis elements in basis. Returns null (0) if multiple*basisElement
    is not regular top reducible -- this indicates a singular
    reduction. */
  virtual Poly* regularReduce(
    const_monomial sig,
    const_monomial multiple,
    size_t basisElement,
    const SigPolyBasis& basis) = 0;

  /** Sets how many bytes of memory to increase the memory use by
    at a time - if such a thing is appropriate for the reducer. */
  virtual void setMemoryQuantum(size_t quantum) = 0;

  // ***** Kinds of reducers and creating a Reducer 

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
    Reducer_Geobucket_Hashed_Packed,

    Reducer_F4_Old,
    Reducer_F4_New
  };

  static std::unique_ptr<Reducer> makeReducer
    (ReducerType t, const PolyRing& ring);

  static std::unique_ptr<Reducer> makeReducerNullOnUnknown
    (ReducerType t, const PolyRing& ring);

  static ReducerType reducerType(int typ);
  static void displayReducerTypes(std::ostream& o);


  // ***** Obtaining statistics about the reduction process

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

  Stats sigStats() const {return mSigStats;}
  Stats classicStats() const {return mClassicStats;}


  // ***** Miscellaneous

  virtual std::string description() const = 0;
  virtual size_t getMemoryUse() const = 0;
  virtual void dump() const;

protected:
  Reducer();

  size_t stats_maxsize;
  size_t stats_maxsize_live;
  unsigned long long stats_n_inserts;
  unsigned long long stats_n_compares;

  Stats mSigStats;
  Stats mClassicStats;
};

MATHICGB_NAMESPACE_END
#endif
