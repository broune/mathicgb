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
#include <fstream>
#include <iostream>

MATHICGB_NAMESPACE_BEGIN

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

  mSPairGroupSize("sPairGroupSize",
    "Specifies how many S-pair to reduce at one time. A value of 0 "
    "indicates to use an appropriate default.",
    0),

 mMinMatrixToStore("storeMatrices",
   "If using a matrix-based reducer, store the matrices that are generated in "
   "files named X-1.mat, X-2.mat and so on where X is the project name. Only "
   "matrices with at least as many entries as the parameter are stored. "
   "A value of 0 indicates not to store any matrices.",
   0),

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
  auto basis = MathicIO<>().readBasis(ring, false, in);

  // run algorithm
  const auto reducerType = Reducer::reducerType(mGBParams.mReducer.value());
  std::unique_ptr<Reducer> reducer;
  if (
    reducerType != Reducer::Reducer_F4_Old &&
    reducerType != Reducer::Reducer_F4_New
  ) {
    reducer = Reducer::makeReducer(reducerType, ring);
  } else {
    const auto type = reducerType == Reducer::Reducer_F4_Old ?
      F4Reducer::OldType : F4Reducer::NewType;
    auto f4Reducer = make_unique<F4Reducer>(ring, type);
    if (mMinMatrixToStore.value() > 0)
      f4Reducer->writeMatricesTo(projectName, mMinMatrixToStore);
    reducer = std::move(f4Reducer);
  }

  ClassicGBAlg alg(
    basis,
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

  if (mGBParams.mOutputResult.value())
    {
      // output Groebner basis into .gb file.

      // The stats information is displayed to cout (above),
      // so we disable it here.
      // std::ofstream statsOut((projectName + ".stats").c_str());
      // alg.printStats(statsOut);

      std::string basisFileName = projectName + ".gb";
      FILE* basisOut = std::fopen(basisFileName.c_str(), "w");
      output(basisOut, alg.basis());
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
}

MATHICGB_NAMESPACE_END
