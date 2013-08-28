// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_GB_COMMON_PARAMS_GUARD
#define MATHICGB_GB_COMMON_PARAMS_GUARD

#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

class GBCommonParams {
public:
  GBCommonParams();

  void pushBackParameters(std::vector<mathic::CliParameter*>& parameters);
  void perform();

  mathic::BoolParameter mPreferSparseReducers;
  mathic::BoolParameter mOutputResult;
  mathic::IntegerParameter mSPairQueue;
  mathic::IntegerParameter mBreakAfter;
  mathic::IntegerParameter mPrintInterval;
  mathic::IntegerParameter mMonomialTable;
  mathic::IntegerParameter mMonoLookup;
  mathic::IntegerParameter mReducer;
  mathic::IntegerParameter mMemoryQuantum;
};

MATHICGB_NAMESPACE_END
#endif
