#include "mathicgb/stdinc.h"
#include "CommonParams.hpp"

CommonParams::CommonParams():
  mInputFile("inputFile",
    "The file to read input from.",
    ""),

  mTracingLevel("tracingLevel",
    "How much information to print out about what the program does. No "
    "information is shown if the value is zero. Higher values "
    "result in more information.",
    0),

  mThreadCount("threadCount",
    "Specifies how many threads to use at a time.",
    1)
{
}

void CommonParams::directOptions(
  std::vector<std::string> tokens,
  mathic::CliParser& parser
) {
  if (tokens.size() == 1)
    mInputFile.processArgument(tokens.back());
  if (tokens.size() > 1)
    mathic::reportError("Too many direct options.");
}

void CommonParams::pushBackParameters(
  std::vector<mathic::CliParameter*>& parameters
) {
  parameters.push_back(&mInputFile);
  parameters.push_back(&mTracingLevel);
  parameters.push_back(&mThreadCount);
}

void CommonParams::perform() {
  tracingLevel = mTracingLevel.value();

  // delete the old init object first to make the new one take control.
  mTbbInit.reset();
  std::unique_ptr<tbb::task_scheduler_init> mTbbInit;

  mTbbInit = make_unique<tbb::task_scheduler_init>(
    mThreadCount.value() == 0 ?
      tbb::task_scheduler_init::automatic :
      mThreadCount.value()
  );
}

void CommonParams::registerFileNameExtension(std::string extension) {
  MATHICGB_ASSERT(!extension.empty());
  mExtensions.push_back(std::move(extension));
}

std::string CommonParams::inputFileNameStem() {
  const auto& str = mInputFile.value();
  const auto toStrip = inputFileNameExtension();
  MATHICGB_ASSERT
    (toStrip.size() < str.size() || (toStrip.empty() && str.empty()));
  return str.substr(0, str.size() - toStrip.size());
}

std::string CommonParams::inputFileNameExtension() {
  const auto& str = mInputFile.value();
  const auto end = mExtensions.end();
  for (auto it = mExtensions.begin(); it != end; ++it) {
    if (
      str.size() >= it->size() &&
      str.substr(str.size() - it->size(), it->size()) == *it
    )
      return *it;
  }
  return std::string();
}
