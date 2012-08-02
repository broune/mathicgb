// Copyright 2011 Michael E. Stillman

#include <iostream>
#include "stdinc.h"
#include "Reducer.hpp"
#include "PolyHeap.hpp"
#include "PolyGeoBucket.hpp"
#include "PolyReducer.hpp"
#include "PolyHashReducer.hpp"
#include "BjarkeGeobucket2.hpp"
#include "GroebnerBasis.hpp"
#include "TournamentReducer.hpp"
#include "HashTourReducer.hpp"
#include "ReducerPack.hpp"
#include "ReducerPackDedup.hpp"
#include "ReducerNoDedup.hpp"
#include "ReducerDedup.hpp"
#include "ReducerHash.hpp"
#include "ReducerHashPack.hpp"
#include <algorithm>

extern int tracingLevel;

std::auto_ptr<Reducer> Reducer::makeReducer(
  ReducerType type,
  PolyRing const& ring
) {
  std::auto_ptr<Reducer> reducer = makeReducerNullOnUnknown(type, ring);
  if (reducer.get() == 0) {
    std::ostringstream error;
    error << "Unknown or unimplemented reducer type " << type << ".\n";
    throw std::runtime_error(error.str());
  }
  return reducer;
}

std::auto_ptr<Reducer> Reducer::makeReducerNullOnUnknown(
  ReducerType type,
  PolyRing const& ring
) {
  switch (type) {
  case Reducer_PolyHeap:
    return std::auto_ptr<Reducer>(new PolyHeap(&ring));
  case Reducer_PolyGeoBucket:
    return std::auto_ptr<Reducer>(new PolyGeoBucket(&ring));
  case Reducer_Poly:
    return std::auto_ptr<Reducer>(new PolyReducer(&ring));
  case Reducer_PolyHash:
    return std::auto_ptr<Reducer>(new PolyHashReducer(&ring));
  case Reducer_BjarkeGeo:
    return std::auto_ptr<Reducer>(new BjarkeGeobucket2(&ring));
  case Reducer_TournamentTree:
    return std::auto_ptr<Reducer>(new TournamentReducer(ring));
    //return std::auto_ptr<Reducer>
      //(new ReducerPack<mic::TourTree>(ring));
  case Reducer_HashTourTree:
    return std::auto_ptr<Reducer>(new HashTourReducer(ring));


  case Reducer_TourTree_NoDedup:
    return std::auto_ptr<Reducer>(new ReducerNoDedup<mic::TourTree>(ring));
  case Reducer_TourTree_Dedup:
    return std::auto_ptr<Reducer>(new ReducerDedup<mic::TourTree>(ring));
  case Reducer_TourTree_Hashed:
    return std::auto_ptr<Reducer>(new ReducerHash<mic::TourTree>(ring));
    //break;
  case Reducer_TourTree_NoDedup_Packed:
    return std::auto_ptr<Reducer>(new ReducerPack<mic::TourTree>(ring));
  case Reducer_TourTree_Dedup_Packed:
    return std::auto_ptr<Reducer>(new ReducerPackDedup<mic::TourTree>(ring));
  case Reducer_TourTree_Hashed_Packed:
    return std::auto_ptr<Reducer>(new ReducerHashPack<mic::TourTree>(ring));

  case Reducer_Heap_NoDedup:
    return std::auto_ptr<Reducer>(new ReducerNoDedup<mic::Heap>(ring));
  case Reducer_Heap_Dedup:
    return std::auto_ptr<Reducer>(new ReducerDedup<mic::Heap>(ring));
  case Reducer_Heap_Hashed:
    return std::auto_ptr<Reducer>(new ReducerHash<mic::Heap>(ring));
    //break;
  case Reducer_Heap_NoDedup_Packed:
    return std::auto_ptr<Reducer>(new ReducerPack<mic::Heap>(ring));
  case Reducer_Heap_Dedup_Packed:
    return std::auto_ptr<Reducer>(new ReducerPackDedup<mic::Heap>(ring));
  case Reducer_Heap_Hashed_Packed:
    return std::auto_ptr<Reducer>(new ReducerHashPack<mic::Heap>(ring));

  case Reducer_Geobucket_NoDedup:
    return std::auto_ptr<Reducer>(new ReducerNoDedup<mic::Geobucket>(ring));
  case Reducer_Geobucket_Dedup:
    return std::auto_ptr<Reducer>(new ReducerDedup<mic::Geobucket>(ring));
  case Reducer_Geobucket_Hashed:
    return std::auto_ptr<Reducer>(new ReducerHash<mic::Geobucket>(ring));
  case Reducer_Geobucket_NoDedup_Packed:
    return std::auto_ptr<Reducer>(new ReducerPack<mic::Geobucket>(ring));
  case Reducer_Geobucket_Dedup_Packed:
    return std::auto_ptr<Reducer>(new ReducerPackDedup<mic::Geobucket>(ring));
  case Reducer_Geobucket_Hashed_Packed:
    return std::auto_ptr<Reducer>(new ReducerHashPack<mic::Geobucket>(ring));

  default:
    break;
  };
  return std::auto_ptr<Reducer>();
}

Reducer::ReducerType Reducer::reducerType(int typ)
{
  switch (typ) {
  case 0: return Reducer_PolyHeap;
  case 1: return Reducer_PolyGeoBucket;
  case 2: return Reducer_Poly;
  case 3: return Reducer_PolyHash;
  case 4: return Reducer_BjarkeGeo;
  case 5: return Reducer_TournamentTree;
  case 6: return Reducer_HashTourTree;

  case 7: return Reducer_TourTree_NoDedup;
  case 8: return Reducer_TourTree_Dedup;
  case 9: return Reducer_TourTree_Hashed;
  case 10: return Reducer_TourTree_NoDedup_Packed;
  case 11: return Reducer_TourTree_Dedup_Packed;
  case 12: return Reducer_TourTree_Hashed_Packed;

  case 13: return Reducer_Heap_NoDedup;
  case 14: return Reducer_Heap_Dedup;
  case 15: return Reducer_Heap_Hashed;
  case 16: return Reducer_Heap_NoDedup_Packed;
  case 17: return Reducer_Heap_Dedup_Packed;
  case 18: return Reducer_Heap_Hashed_Packed;

  case 19: return Reducer_Geobucket_NoDedup;
  case 20: return Reducer_Geobucket_Dedup;
  case 21: return Reducer_Geobucket_Hashed;
  case 22: return Reducer_Geobucket_NoDedup_Packed;
  case 23: return Reducer_Geobucket_Dedup_Packed;
  case 24: return Reducer_Geobucket_Hashed_Packed;

  default: return Reducer_PolyHeap;
  }
}

void Reducer::displayReducerTypes(std::ostream &o)
{
  o << "Reducer types:" << std::endl;
  o << "   0   PolyHeap" << std::endl;
  o << "   1   PolyGeoBucket" << std::endl;
  o << "   2   Poly" << std::endl;
  o << "   3   PolyHash" << std::endl;
  o << "   4   BjarkeGeo2" << std::endl;
  o << "   5   Tournament tree" << std::endl;
  o << "   6   Hashed Tournament tree" << std::endl;

  o << "   7   TournamentTree.NoDedup" << std::endl;
  o << "   8   TournamentTree.Dedup" << std::endl;
  o << "   9   TournamentTree.Hashed" << std::endl;
  o << "  10   TournamentTree.NoDedup.Packed" << std::endl;
  o << "  11   TournamentTree.Dedup.Packed" << std::endl;
  o << "  12   TournamentTree.Hashed.Packed" << std::endl;

  o << "  13   Heap.NoDedup" << std::endl; 
  o << "  14   Heap.Dedup" << std::endl;
  o << "  15   Heap.Hashed" << std::endl;
  o << "  16   Heap.NoDedup.Packed" << std::endl;
  o << "  17   Heap.Dedup.Packed" << std::endl;
  o << "  18   Heap.Hashed.Packed" << std::endl;

  o << "  19   Geobucket.NoDedup" << std::endl;
  o << "  20   Geobucket.Dedup" << std::endl;
  o << "  21   Geobucket.Hashed" << std::endl;
  o << "  22   Geobucket.NoDedup.Packed" << std::endl;
  o << "  23   Geobucket.Dedup.Packed" << std::endl;
  o << "  24   Geobucket.Hashed.Packed" << std::endl;
}

///////////////////////////
// ReducerPack //////////
///////////////////////////

///////////////////////////
// Reducer NoDedup/Dedup //
///////////////////////////


///////////////////////
// Reducer ////////////
/////////////////////// 

Reducer::Stats::Stats():
  reductions(0),
  singularReductions(0),
  zeroReductions(0),
  steps(0),
  maxSteps(0) {}

Reducer::Reducer():
  stats_maxsize(0),
  stats_maxsize_live(0),
  stats_n_inserts(0),
  stats_n_compares(0) {}

void Reducer::reset()
{
  mArena.freeAllAllocs();
  resetReducer();
}

size_t Reducer::getMemoryUse() const {
  return mArena.getMemoryUse();
}

Poly* Reducer::regularReduce(
  const_monomial sig,
  const_monomial multiple,
  size_t basisElement,
  const GroebnerBasis& basis)
{
  const PolyRing& ring = basis.ring();
  ++mSigStats.reductions;

  monomial tproduct = ring.allocMonomial(mArena);
  monomial u = ring.allocMonomial(mArena);
  ring.monomialMult(multiple, basis.getLeadMonomial(basisElement), tproduct);

  size_t reducer = basis.regularReducer(sig, tproduct);
  if (reducer == static_cast<size_t>(-1)) {
    ++mSigStats.singularReductions;
    mArena.freeAllAllocs();
    return 0; // singular reduction: no regular top reduction possible
  }

  ring.monomialDivide(tproduct, basis.getLeadMonomial(reducer), u);

  coefficient coef;
  ring.coefficientSet(coef, 1);
  insertTail(const_term(coef, multiple), &basis.poly(basisElement));

  ASSERT(ring.coefficientIsOne(basis.getLeadCoefficient(reducer)));
  ring.coefficientFromInt(coef, -1);
  insertTail(const_term(coef, u), &basis.poly(reducer));
  basis.basis().usedAsReducer(reducer);

  Poly* result = new Poly(&ring);

  unsigned long long steps = 2; // number of steps in this reduction
  for (const_term v; findLeadTerm(v); ++steps) {
    ASSERT(v.coeff != 0);
    reducer = basis.regularReducer(sig, v.monom);
    if (reducer == static_cast<size_t>(-1)) { // no reducer found
      result->appendTerm(v.coeff, v.monom);
      removeLeadTerm();
    } else { // reduce by reducer
      basis.basis().usedAsReducer(reducer);
      monomial mon = ring.allocMonomial(mArena);
      ring.monomialDivide(v.monom, basis.getLeadMonomial(reducer), mon);
      ring.coefficientDivide(v.coeff, basis.getLeadCoefficient(reducer), coef);
      ring.coefficientNegateTo(coef);
      removeLeadTerm();
      insertTail(const_term(coef, mon), &basis.poly(reducer));
    }
  }
  result->makeMonic();

  mSigStats.steps += steps;
  mSigStats.maxSteps = std::max(mSigStats.maxSteps, steps);
  if (result->isZero())
    ++mSigStats.zeroReductions;

  reset();
  return result;
}

std::auto_ptr<Poly> Reducer::classicReduce(const Poly& poly, const PolyBasis& basis) {
  monomial identity = basis.ring().allocMonomial(mArena);
  basis.ring().monomialSetIdentity(identity);
  insert(identity, &poly);

  return classicReduce(basis);
}

std::auto_ptr<Poly> Reducer::classicTailReduce(const Poly& poly, const PolyBasis& basis) {
  ASSERT(&poly.ring() == &basis.ring());
  ASSERT(!poly.isZero());
  term identity;
  identity.monom = basis.ring().allocMonomial(mArena);
  basis.ring().monomialSetIdentity(identity.monom);
  basis.ring().coefficientSetOne(identity.coeff);
  insertTail(identity, &poly);

  std::auto_ptr<Poly> result(new Poly(&basis.ring()));
  result->appendTerm(poly.getLeadCoefficient(), poly.getLeadMonomial());

  return classicReduce(result, basis);
}

std::auto_ptr<Poly> Reducer::classicReduceSPoly(
  const Poly& a,
  const Poly& b,
  const PolyBasis& basis
) {
  const PolyRing& ring = basis.ring();

  monomial lcm = ring.allocMonomial();
  ring.monomialLeastCommonMultiple
    (a.getLeadMonomial(), b.getLeadMonomial(), lcm);

  // insert tail of multiple of a
  monomial multiple1 = ring.allocMonomial();
  ring.monomialDivide(lcm, a.getLeadMonomial(), multiple1);
  coefficient plusOne;
  ring.coefficientSet(plusOne, 1);
  insertTail(const_term(plusOne, multiple1), &a);

  // insert tail of multiple of b
  monomial multiple2 = ring.allocMonomial();
  ring.monomialDivide(lcm, b.getLeadMonomial(), multiple2);
  coefficient minusOne = plusOne;
  ring.coefficientNegateTo(minusOne);
  insertTail(const_term(minusOne, multiple2), &b);

  std::auto_ptr<Poly> reduced = classicReduce(basis);
  ring.freeMonomial(lcm);
  ring.freeMonomial(multiple1);
  ring.freeMonomial(multiple2);
  return reduced;
}

std::auto_ptr<Poly> Reducer::classicReduce
    (std::auto_ptr<Poly> result, const PolyBasis& basis) {
  const PolyRing& ring = basis.ring();
  ASSERT(&result->ring() == &ring);
  ++mClassicStats.reductions;

  if (tracingLevel > 100)
    std::cerr << "Classic reduction begun." << std::endl;

  coefficient coef;
  unsigned long long steps = 1; // number of steps in this reduction
  for (const_term v; findLeadTerm(v); ++steps) {
    if (tracingLevel > 100) {
      std::cerr << "from reducer queue: ";
      basis.ring().monomialDisplay(std::cerr, v.monom);
      std::cerr << std::endl;
    }

    size_t reducer = basis.divisor(v.monom);
    if (reducer == static_cast<size_t>(-1)) { // no reducer found
      ASSERT(result->isZero() ||
        basis.order().signatureCompare(v.monom, result->backMonomial()) == LT);
      result->appendTerm(v.coeff, v.monom);
      removeLeadTerm();
    } else { // reduce by reducer
      basis.usedAsReducer(reducer);
      monomial mon = ring.allocMonomial(mArena);
      ring.monomialDivide(v.monom, basis.leadMonomial(reducer), mon);
      ring.coefficientDivide(v.coeff, basis.leadCoefficient(reducer), coef);
      ring.coefficientNegateTo(coef);
      removeLeadTerm();
      insertTail(const_term(coef, mon), &basis.poly(reducer));

      if (tracingLevel > 100) {
        std::cerr << "Reducing by basis element " << reducer << ": ";
        basis.poly(reducer).display(std::cerr);
        std::cerr << std::endl;
        std::cerr << "multiplied by: " << coef << "  *  ";
        basis.ring().monomialDisplay(std::cerr, mon);
        std::cerr << std::endl;
      }
    }
  }
  result->makeMonic();

  mClassicStats.steps += steps;
  mClassicStats.maxSteps = std::max(mClassicStats.maxSteps, steps);
  if (result->isZero())
    ++mClassicStats.zeroReductions;

  if (tracingLevel > 100)
    std::cerr << "Classic reduction done." << std::endl;

  reset();
  return result;
}

std::auto_ptr<Poly> Reducer::classicReduce(const PolyBasis& basis) {
  std::auto_ptr<Poly> result(new Poly(&basis.ring()));
  return classicReduce(result, basis);
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
