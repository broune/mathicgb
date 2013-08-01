// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_HELP_ACTION_GUARD
#define MATHICGB_HELP_ACTION_GUARD

#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

class HelpAction : public mathic::HelpAction {
public:
  virtual void performAction();
};

MATHICGB_NAMESPACE_END

#endif
