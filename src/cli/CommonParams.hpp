// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_COMMON_PARAMS_GUARD
#define MATHICGB_COMMON_PARAMS_GUARD

#include "mathicgb/mtbb.hpp"
#include <mathic.h>
#include <vector>

MATHICGB_NAMESPACE_BEGIN

class CommonParams {
public:
  CommonParams(size_t minDirectParams, size_t maxDirectParams);

  void directOptions
    (std::vector<std::string> tokens, mathic::CliParser& parser);
    
  void pushBackParameters(std::vector<mathic::CliParameter*>& parameters);

  /// Takes appropriate action depending on the parameters. For example this
  /// will set the number of threads in tbb.
  void perform();

  /// If called with string X, then X will be considered an extension
  /// for a file name instead of part of the file name.
  void registerFileNameExtension(std::string extensions);

  /// Returns the number of direct parameters/input files.
  size_t inputFileCount() const;

  /// Returns the file name at offset i, if any.
  std::string inputFileName(size_t i);

  /// Returns the stem of the input file name at offset i, with any registered
  /// extensions stripped off.
  std::string inputFileNameStem(size_t i);

  /// Returns the registered extension of the input file name at offset i,
  /// if any.
  std::string inputFileNameExtension(size_t i);

private:
  mathic::IntegerParameter mTracingLevel;
  mathic::IntegerParameter mThreadCount;
  mathic::StringParameter mLogs;

  std::vector<std::string> mExtensions; /// to recognize file type

  /// to set thread count
  std::unique_ptr<mgb::mtbb::task_scheduler_init> mTbbInit;
  std::size_t mMinDirectParams;
  std::size_t mMaxDirectParams;
  std::vector<std::string> mDirectParameters;
};

MATHICGB_NAMESPACE_END

#endif
