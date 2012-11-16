#ifndef MATHICGB_COMMON_PARAMS_GUARD
#define MATHICGB_COMMON_PARAMS_GUARD

#include <mathic.h>
#include <tbb/tbb.h>
#include <vector>

class CommonParams {
public:
  CommonParams();

  void directOptions
    (std::vector<std::string> tokens, mathic::CliParser& parser);
    
  void pushBackParameters(std::vector<mathic::CliParameter*>& parameters);

  /// Takes appropriate action depending on the parameters. For example this
  /// will set the number of threads in tbb.
  void perform();

  /// If called with string X, then X will be considered an extension
  /// for a file name instead of part of the file name.
  void registerFileNameExtension(std::string extensions);

  /// Returns the stem of the input file name, with any registered extensions
  /// stripped off.
  std::string inputFileNameStem();

  /// Returns the registered extension of the input file name, if any.
  std::string inputFileNameExtension();

  mathic::StringParameter mInputFile;
  mathic::IntegerParameter mTracingLevel;
  mathic::IntegerParameter mThreadCount;

private:
  std::vector<std::string> mExtensions; // to recognize file type
  std::unique_ptr<tbb::task_scheduler_init> mTbbInit; // to set thread count
};

#endif
