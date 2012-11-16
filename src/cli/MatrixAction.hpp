#ifndef MATHICGB_MATRIX_ACTION_GUARD
#define MATHICGB_MATRIX_ACTION_GUARD

#include "CommonParams.hpp"
#include <mathic.h>

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

#endif
