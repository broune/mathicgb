// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_BJARKE_GEOBUCKET2_GUARD
#define MATHICGB_BJARKE_GEOBUCKET2_GUARD

MATHICGB_NAMESPACE_BEGIN

class TypicalReducer;
class PolyRing;

std::unique_ptr<TypicalReducer> makeBjarkeGeobucket2(const PolyRing& ring);

MATHICGB_NAMESPACE_END
#endif
