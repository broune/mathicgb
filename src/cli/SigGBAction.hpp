#ifndef MATHICGB_SIG_G_B_ACTION_GUARD
#define MATHICGB_SIG_G_B_ACTION_GUARD

#include "GBCommonParams.hpp"
#include "CommonParams.hpp"
#include <mathic.h>

class SigGBAction : public mic::Action {
public:
  SigGBAction();

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
  mic::BoolParameter mUseSingularCriterionEarly;
  mic::BoolParameter mPostponeKoszul;
  mic::BoolParameter mUseBaseDivisors;
  mic::IntegerParameter mModuleOrder;
};

#endif
