// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_MATRIX_ACTION_GUARD
#define MATHICGB_MATRIX_ACTION_GUARD

#include "CommonParams.hpp"
#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

/// Performs computations on matrices.
class MatrixAction : public mathic::Action {
public:
  MatrixAction();

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
};

MATHICGB_NAMESPACE_END
#endif
