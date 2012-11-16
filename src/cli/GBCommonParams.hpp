#ifndef MATHICGB_GB_COMMON_PARAMS_GUARD
#define MATHICGB_GB_COMMON_PARAMS_GUARD

#include <mathic.h>

class GBCommonParams {
public:
  GBCommonParams();

  void pushBackParameters(std::vector<mathic::CliParameter*>& parameters);
  void perform();

  mathic::BoolParameter mPreferSparseReducers;
  mathic::IntegerParameter mSPairQueue;
  mathic::IntegerParameter mBreakAfter;
  mathic::IntegerParameter mPrintInterval;
  mathic::IntegerParameter mMonomialTable;
  mathic::IntegerParameter mDivisorLookup;
  mathic::IntegerParameter mReducer;
  mathic::IntegerParameter mMemoryQuantum;
};

#endif
