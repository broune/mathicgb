// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "mathicgb/stdinc.h"
#include "GBAction.hpp"

#include "mathicgb/ClassicGBAlg.hpp"
#include "mathicgb/Basis.hpp"
#include "mathicgb/io-util.hpp"
#include "mathicgb/F4Reducer.hpp"
#include "mathicgb/Scanner.hpp"
#include "mathicgb/MathicIO.hpp"
#include "mathicgb/Reducer.hpp"
#include <fstream>
#include <iostream>

MATHICGB_NAMESPACE_BEGIN

GBAction::GBAction():
  mAutoTailReduce(
    "autoTailReduce",
    "Reduce the non-leading terms of all polynomials whenever an element "
    "is inserted into the basis. Only relevant to the "
    "classic Buchberger algorithm.",
    false),

  mAutoTopReduce(
    "autoTopReduce",
    "Reduce any basis element whose lead term becomes reducible "
    "by a different basis element. Only relevant to the "
    "classic Buchberger algorithm.",
    true),

  mSPairGroupSize(
    "sPairGroupSize",
    "Specifies how many S-pair to reduce at one time. A value of 0 "
    "indicates to use an appropriate default.",
    0),

  mMinMatrixToStore(
    "storeMatrices",
    "If using a matrix-based reducer, store the matrices that are generated in "
    "files named X-1.mat, X-2.mat and so on where X is the project name. Only "
    "matrices with at least as many entries as the parameter are stored. "
    "A value of 0 indicates not to store any matrices.",
    0),

  mModule(
    "module",
    "The input is a basis of a submodule over the polynomial ring instead of "
    "an ideal in the polynomial ring. This option is experimental.",
    false
  ),

   mParams(1, 1)
{}

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
  const std::string inputBasisFile = projectName + ".ideal";
  std::ifstream inputFile(inputBasisFile.c_str());
  if (inputFile.fail())
    mic::reportError("Could not read input file \"" + inputBasisFile + '\n');

  Scanner in(inputFile);
  auto p = MathicIO<>().readRing(true, in);
  auto& ring = *p.first;
  auto basis = MathicIO<>().readBasis(ring, mModule.value(), in);

  // run algorithm
  const auto reducerType = Reducer::reducerType(mGBParams.mReducer.value());
  std::unique_ptr<Reducer> reducer;
  if (
    reducerType != Reducer::Reducer_F4_Old &&
    reducerType != Reducer::Reducer_F4_New
  ) {
    reducer = Reducer::makeReducer(reducerType, ring);
  } else {
    auto f4Reducer = makeF4Reducer(
      ring,
      reducerType == Reducer::Reducer_F4_Old,
      mMinMatrixToStore.value() > 0 ? projectName : "",
      mMinMatrixToStore
    );     
    reducer = std::move(f4Reducer);
  }

  ClassicGBAlgParams params;
  params.reducer = reducer.get();
  params.monoLookupType = mGBParams.mMonoLookup.value();
  params.preferSparseReducers = mGBParams.mPreferSparseReducers.value();
  params.sPairQueueType = mGBParams.mSPairQueue.value();
  params.breakAfter = mGBParams.mBreakAfter.value();
  params.printInterval = mGBParams.mPrintInterval.value();
  params.sPairGroupSize = mSPairGroupSize.value();
  params.reducerMemoryQuantum = mGBParams.mMemoryQuantum.value();
  params.useAutoTopReduction = mAutoTopReduce.value();
  params.useAutoTailReduction = mAutoTailReduce.value();
  params.callback = nullptr;

  const auto gb = mModule.value() ?
    computeModuleGBClassicAlg(std::move(basis), params) :
    computeGBClassicAlg(std::move(basis), params);

  if (mGBParams.mOutputResult.value()) {
    std::ofstream out(projectName + ".gb");
    MathicIO<>().writeBasis(gb, mModule.value(), out);
  }
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
  parameters.push_back(&mModule);
}

MATHICGB_NAMESPACE_END
  