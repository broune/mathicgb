#include "mathicgb/stdinc.h"
#include "GBAction.hpp"

#include "mathicgb/BuchbergerAlg.hpp"
#include "mathicgb/Ideal.hpp"
#include "mathicgb/io-util.hpp"
#include "mathicgb/F4Reducer.hpp"
#include <fstream>
#include <iostream>

GBAction::GBAction():
  mAutoTailReduce("autoTailReduce",
    "Reduce the non-leading terms of all polynomials whenever an element "
    "is inserted into the basis. Only relevant to the "
    "classic Buchberger algorithm.",
    false),

  mAutoTopReduce("autoTopReduce",
    "Reduce any basis element whose lead term becomes reducible "
    "by a different basis element. Only relevant to the "
    "classic Buchberger algorithm.",
    true),

  //mTermOrder("termOrder",
  //  "The term order to compute a Grobner basis with respect to. This is "
  //  "currently actually a free module term order, but that should be changed.",
  //  4),

  mSPairGroupSize("sPairGroupSize",
    "Specifies how many S-pair to reduce at one time. A value of 0 "
    "indicates not to group S-pairs together. Only currently relevant "
    "for the classic Buchberger algorithm.",
    0),

 mMinMatrixToStore("storeMatrices",
   "If using a matrix-based reducer, store the matrices that are generated in "
   "files named X-1.mat, X-2.mat and so on where X is the project name. Only "
   "matrices with at least as many entries as the parameter are stored. "
   "A value of 0 indicates not to store any matrices.",
   0),

   mParams(1, 1)
{
  std::ostringstream orderOut;
  FreeModuleOrder::displayOrderTypes(orderOut);
  //mTermOrder.appendToDescription(orderOut.str());
}

void GBAction::directOptions(
  std::vector<std::string> tokens,
  mic::CliParser& parser
) {
  mParams.directOptions(tokens, parser);
}

void GBAction::performAction() {
  mParams.perform();
  mGBParams.perform();
  const std::string projectName = mParams.inputFileNameStem(0);

  // read input
  std::unique_ptr<Ideal> ideal;
  {
    const std::string inputIdealFile = projectName + ".ideal";
    std::ifstream inputFile(inputIdealFile.c_str());
    if (inputFile.fail())
      mic::reportError("Could not read input file \"" + inputIdealFile + '\n');
    ideal = Ideal::parse(inputFile);
  }
  std::unique_ptr<PolyRing const> ring(&(ideal->ring()));

  // run algorithm
  const auto reducerType = Reducer::reducerType(mGBParams.mReducer.value());
  std::unique_ptr<Reducer> reducer;
  if (reducerType != Reducer::Reducer_F4)
    reducer = Reducer::makeReducer(reducerType, *ring);
  else {
    auto f4Reducer = make_unique<F4Reducer>(*ring);
    if (mMinMatrixToStore.value() > 0)
      f4Reducer->writeMatricesTo(projectName, mMinMatrixToStore);
    reducer = std::move(f4Reducer);
  }

  BuchbergerAlg alg(
    *ideal,
    4 /*mModuleOrder.value()*/ , // todo: alg should not take a *module* order
    *reducer,
    mGBParams.mDivisorLookup.value(),
    mGBParams.mPreferSparseReducers.value(),
    mGBParams.mSPairQueue.value());
  alg.setBreakAfter(mGBParams.mBreakAfter.value());
  alg.setPrintInterval(mGBParams.mPrintInterval.value());
  alg.setSPairGroupSize(mSPairGroupSize.value());
  alg.setReducerMemoryQuantum(mGBParams.mMemoryQuantum.value());
  alg.setUseAutoTopReduction(mAutoTopReduce.value());
  alg.setUseAutoTailReduction(mAutoTailReduce.value());

  alg.computeGrobnerBasis();
  alg.printStats(std::cerr);

  /*
  std::ofstream statsOut((mProjectName.value() + ".stats").c_str());
  alg.printStats(statsOut);
  {
    std::string basisFileName = mProjectName.value() + ".gb";
    FILE* basisOut = std::fopen(basisFileName.c_str(), "w");
    output(basisOut, alg.basis());
  }
  */
}

const char* GBAction::staticName() {
  return "gb";
}

const char* GBAction::name() const {
  return staticName();
}

const char* GBAction::description() const {
  return "Compute a Grobner basis. "
    "The project name is an optional direct parameter.";
}

const char* GBAction::shortDescription() const {
  return "Compute a Grobner basis.";
}

void GBAction::pushBackParameters(
  std::vector<mic::CliParameter*>& parameters
) {
  mParams.pushBackParameters(parameters);
  mGBParams.pushBackParameters(parameters);
  parameters.push_back(&mAutoTailReduce);
  parameters.push_back(&mAutoTopReduce);
  parameters.push_back(&mSPairGroupSize);
  parameters.push_back(&mMinMatrixToStore);
}
