#include "stdinc.h"
#include "Reducer.hpp"

#include "PolyHeap.hpp"
#include "PolyGeoBucket.hpp"
#include "PolyReducer.hpp"
#include "PolyHashReducer.hpp"
#include "BjarkeGeobucket2.hpp"
#include "TournamentReducer.hpp"
#include "HashTourReducer.hpp"
#include "ReducerPack.hpp"
#include "ReducerPackDedup.hpp"
#include "ReducerNoDedup.hpp"
#include "ReducerDedup.hpp"
#include "ReducerHash.hpp"
#include "ReducerHashPack.hpp"
#include "F4Reducer.hpp"

#include "SigPolyBasis.hpp"
#include <iostream>
#include <algorithm>

Reducer::Reducer():
  stats_maxsize(0),
  stats_maxsize_live(0),
  stats_n_inserts(0),
  stats_n_compares(0) {}

Reducer::~Reducer() {
}

std::unique_ptr<Reducer> Reducer::makeReducer(
  ReducerType type,
  PolyRing const& ring
) {
  std::unique_ptr<Reducer> reducer =
    makeReducerNullOnUnknown(type, ring);
  if (reducer.get() == 0) {
    std::ostringstream error;
    error << "Unknown or unimplemented reducer type " << type << ".\n";
    throw std::runtime_error(error.str());
  }
  return reducer;
}

std::unique_ptr<Reducer> Reducer::makeReducerNullOnUnknown(
  ReducerType type,
  PolyRing const& ring
) {
  switch (type) {
  case Reducer_PolyHeap:
    return std::unique_ptr<Reducer>(new PolyHeap(&ring));
  case Reducer_PolyGeoBucket:
    return std::unique_ptr<Reducer>(new PolyGeoBucket(&ring));
  case Reducer_Poly:
    return std::unique_ptr<Reducer>(new PolyReducer(&ring));
  case Reducer_PolyHash:
    return std::unique_ptr<Reducer>(new PolyHashReducer(&ring));
  case Reducer_BjarkeGeo:
    return std::unique_ptr<Reducer>(new BjarkeGeobucket2(&ring));
  case Reducer_TournamentTree:
    return std::unique_ptr<Reducer>(new TournamentReducer(ring));
  case Reducer_HashTourTree:
    return std::unique_ptr<Reducer>(new HashTourReducer(ring));

  case Reducer_TourTree_NoDedup:
    return std::unique_ptr<Reducer>(new ReducerNoDedup<mic::TourTree>(ring));
  case Reducer_TourTree_Dedup:
    return std::unique_ptr<Reducer>(new ReducerDedup<mic::TourTree>(ring));
  case Reducer_TourTree_Hashed:
    return std::unique_ptr<Reducer>(new ReducerHash<mic::TourTree>(ring));
  case Reducer_TourTree_NoDedup_Packed:
    return std::unique_ptr<Reducer>(new ReducerPack<mic::TourTree>(ring));
  case Reducer_TourTree_Dedup_Packed:
    return std::unique_ptr<Reducer>(new ReducerPackDedup<mic::TourTree>(ring));
  case Reducer_TourTree_Hashed_Packed:
    return std::unique_ptr<Reducer>(new ReducerHashPack<mic::TourTree>(ring));

  case Reducer_Heap_NoDedup:
    return std::unique_ptr<Reducer>(new ReducerNoDedup<mic::Heap>(ring));
  case Reducer_Heap_Dedup:
    return std::unique_ptr<Reducer>(new ReducerDedup<mic::Heap>(ring));
  case Reducer_Heap_Hashed:
    return std::unique_ptr<Reducer>(new ReducerHash<mic::Heap>(ring));
  case Reducer_Heap_NoDedup_Packed:
    return std::unique_ptr<Reducer>(new ReducerPack<mic::Heap>(ring));
  case Reducer_Heap_Dedup_Packed:
    return std::unique_ptr<Reducer>(new ReducerPackDedup<mic::Heap>(ring));
  case Reducer_Heap_Hashed_Packed:
    return std::unique_ptr<Reducer>(new ReducerHashPack<mic::Heap>(ring));

  case Reducer_Geobucket_NoDedup:
    return std::unique_ptr<Reducer>(new ReducerNoDedup<mic::Geobucket>(ring));
  case Reducer_Geobucket_Dedup:
    return std::unique_ptr<Reducer>(new ReducerDedup<mic::Geobucket>(ring));
  case Reducer_Geobucket_Hashed:
    return std::unique_ptr<Reducer>(new ReducerHash<mic::Geobucket>(ring));
  case Reducer_Geobucket_NoDedup_Packed:
    return std::unique_ptr<Reducer>(new ReducerPack<mic::Geobucket>(ring));
  case Reducer_Geobucket_Dedup_Packed:
    return std::unique_ptr<Reducer>(new ReducerPackDedup<mic::Geobucket>(ring));
  case Reducer_Geobucket_Hashed_Packed:
    return std::unique_ptr<Reducer>(new ReducerHashPack<mic::Geobucket>(ring));

  case Reducer_F4_Old:
    return make_unique<F4Reducer>(ring, F4Reducer::OldType);
  case Reducer_F4_New:
    return make_unique<F4Reducer>(ring, F4Reducer::NewType);

  default:
    break;
  };
  return std::unique_ptr<Reducer>();
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

  case 25: return Reducer_F4_Old;
  case 26: return Reducer_F4_New;

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

  o << "  25   F4 reducer, old" << std::endl;
  o << "  26   F4 reducer, new" << std::endl;
}

// Todo: can't this be machine generated?
Reducer::Stats::Stats():
  reductions(0),
  singularReductions(0),
  zeroReductions(0),
  steps(0),
  maxSteps(0) {}

void Reducer::dump() const {
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
