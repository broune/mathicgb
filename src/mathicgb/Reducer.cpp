// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "Reducer.hpp"

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

MATHICGB_NAMESPACE_BEGIN

// Calling these dummy methods from the various reducer file headers ensures
// that those translation units are registered as necessary. Without this,
// they will not be linked in, because apparently the linker does not consider
// the presense of global object constructors as a reason to include a
// translation unit into the library. It will simply note that no other
// translation unit has a dependency on anything in that translation unit
// and then not include it in the library. This was wonderful to diagnose.
void dummyLinkerFix() {
  reducerPackDependency();
}

Reducer::~Reducer() {
  // Has to be called somewhere or GCC will eliminate it.
  dummyLinkerFix();
}

/// Vector that stores the registered reducer typers. This has to be a
/// function rather than just a naked object to ensure that the object
/// gets initialized before it is used.
std::vector<Reducer::Registration*>& reducerTypes() {
  static std::vector<Reducer::Registration*> types;
  return types;
}

Reducer::Registration::Registration(
  const char* name, 
  ReducerType id,
  std::unique_ptr<Reducer> (*create)(const PolyRing&)
):
  mName(name),
  mId(id),
  mCreate(create)
{
  reducerTypes().push_back(this);
}

Reducer::Registration r1(
  "TourNoDedup",
  Reducer::Reducer_TourTree_NoDedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerNoDedup<mic::TourTree>>(ring);
  }
);

Reducer::Registration r2(
  "TourDedup",
  Reducer::Reducer_TourTree_Dedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerDedup<mic::TourTree>>(ring);
  }
);
Reducer::Registration r3(
  "TourHash",
  Reducer::Reducer_TourTree_Hashed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHash<mic::TourTree>>(ring);
  }
);
 
 
Reducer::Registration r5(
  "TourDedupPack",
  Reducer::Reducer_TourTree_Dedup_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerPackDedup<mic::TourTree>>(ring);
  }
);
  
Reducer::Registration r6(
  "TourHashPack",
  Reducer::Reducer_TourTree_Hashed_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHashPack<mic::TourTree>>(ring);
  }
);

Reducer::Registration r7(
  "HeapNoDedup",
  Reducer::Reducer_Heap_NoDedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerNoDedup<mic::Heap>>(ring);
  }
);
Reducer::Registration r8(
  "HeapDedup",
  Reducer::Reducer_Heap_Dedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerDedup<mic::Heap>>(ring);
  }
);
 
Reducer::Registration r9(
  "HeapHash",
  Reducer::Reducer_Heap_Hashed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHash<mic::Heap>>(ring);
  }
);
Reducer::Registration r11(
  "HeapDedupPack",
  Reducer::Reducer_Heap_Dedup_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerPackDedup<mic::Heap>>(ring);
  }
);
Reducer::Registration r12(
  "HeapHashPack",
  Reducer::Reducer_Heap_Hashed_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHashPack<mic::Heap>>(ring);
  }
);
Reducer::Registration r13(
  "GeoNoDedup",
  Reducer::Reducer_Geobucket_NoDedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerNoDedup<mic::Geobucket>>(ring);
  }
);
Reducer::Registration r14(
  "GeoDedup",
  Reducer::Reducer_Geobucket_Dedup,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerDedup<mic::Geobucket>>(ring);
  }
);
Reducer::Registration r15(
  "GeoHash",
  Reducer::Reducer_Geobucket_Hashed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHash<mic::Geobucket>>(ring);
  }
);
Reducer::Registration r17(
  "GeoDedupPack",
  Reducer::Reducer_Geobucket_Dedup_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerPackDedup<mic::Geobucket>>(ring);
  }
);
Reducer::Registration r18(
  "GeoHashPack",
  Reducer::Reducer_Geobucket_Hashed_Packed,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<ReducerHashPack<mic::Geobucket>>(ring);
  }
);
Reducer::Registration r19(
  "F4Old",
  Reducer::Reducer_F4_Old,
  [](const PolyRing& ring) -> std::unique_ptr<Reducer> {
    return make_unique<F4Reducer>(ring, F4Reducer::OldType);
  }
);

MATHICGB_REGISTER_REDUCER(
  "F4New",
  Reducer_F4_New,
  (make_unique<F4Reducer>(ring, F4Reducer::NewType))
);

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
  for (const auto& r : reducerTypes()) {
    if (type == r->mId)
      return r->mCreate(ring);
  }
  return nullptr;
}

Reducer::ReducerType Reducer::reducerType(int type)
{
  switch (type) {
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

  default: return Reducer_Geobucket_Hashed;
  }
}

void Reducer::displayReducerTypes(std::ostream& out)
{
  mathic::ColumnPrinter pr;
  auto& id = pr.addColumn(false, "  ");
  auto& desc = pr.addColumn(true, "   ");
  for (const auto& r : reducerTypes()) {
    id << r->mId << '\n';
    desc << r->mName << '\n';
  }
  out << "Reducer types:\n" << pr;
}

MATHICGB_NAMESPACE_END
