#ifndef MATHICGB_G_B_ACTION_GUARD
#define MATHICGB_G_B_ACTION_GUARD

#include "GBCommonParams.hpp"
#include "CommonParams.hpp"
#include <mathic.h>

/// Calculates a classic Grobner basis using Buchberger's algorithm
class GBAction : public mathic::Action {
public:
  GBAction();

  virtual void directOptions(
    std::vector<std::string> tokens,
    mic::CliParser& parser
  );

  virtual void performAction();

  static const char* staticName();

  virtual const char* name() const;
  virtual const char* description() const;
  virtual const char* shortDescription() const;
  
  virtual void pushBackParameters(std::vector<mic::CliParameter*>& parameters);

private:
  CommonParams mParams;
  GBCommonParams mGBParams;
  mathic::BoolParameter mAutoTailReduce;
  mathic::BoolParameter mAutoTopReduce;
  //mic::IntegerParameter mTermOrder;
  mathic::IntegerParameter mSPairGroupSize;
  mathic::IntegerParameter mMinMatrixToStore;
};

#endif
