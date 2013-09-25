// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "mathicgb/stdinc.h"
#include "SigGBAction.hpp"

#include "mathicgb/Basis.hpp"
#include "mathicgb/SignatureGB.hpp"
#include "mathicgb/io-util.hpp"
#include "mathicgb/Scanner.hpp"
#include "mathicgb/MathicIO.hpp"
#include <fstream>
#include <iostream>

MATHICGB_NAMESPACE_BEGIN

SigGBAction::SigGBAction():
  mUseSingularCriterionEarly(
    "earlySingularCriterion",
    "Apply the singular S-pair elimination criterion before queueing "
    "that S-pair. Otherwise, the criterion is only checked just before "
    "the S-pair would otherwise cause a polynomial reduction to occur. ",
    false),

  mPostponeKoszul(
    "postponeKoszul",
    "Postpone the construction of Koszul syzygy signatures.",
    true),

  mUseBaseDivisors(
    "useBaseDivisors",
    "Use high ratio and low ratio base divisors to eliminate "
    "S-spairs quickly based on signature.",
    true),

  mParams(1, 1)
{}

void SigGBAction::directOptions(
  std::vector< std::string> tokens,
  mic::CliParser& parser
) {
  mParams.directOptions(tokens, parser);
}

void SigGBAction::performAction() {
  mParams.perform();
  mGBParams.perform();

  // read input file
  const std::string inputBasisFile = mParams.inputFileNameStem(0) + ".ideal";
  std::ifstream inputFile(inputBasisFile.c_str());
  if (inputFile.fail())
    mic::reportError("Could not read input file \"" + inputBasisFile + '\n');


  Scanner in(inputFile);
  auto p = MathicIO<>().readRing(true, in);
  auto& ring = *p.first;
  auto& processor = p.second;
  auto basis = MathicIO<>().readBasis(ring, false, in);
  if (processor.schreyering())
    processor.setSchreyerMultipliers(basis);

  SignatureGB alg(
    std::move(basis),
    std::move(processor),
    Reducer::reducerType(mGBParams.mReducer.value()),
    mGBParams.mMonoLookup.value(),
    mGBParams.mMonomialTable.value(),
    mPostponeKoszul.value(),
    mUseBaseDivisors.value(),
    mGBParams.mPreferSparseReducers.value(),
    mUseSingularCriterionEarly.value(),
    mGBParams.mSPairQueue.value());
  alg.setBreakAfter(mGBParams.mBreakAfter.value());
  alg.setPrintInterval(mGBParams.mPrintInterval.value());
  alg.computeGrobnerBasis();

  // print statistics
  alg.displayStats(std::cout);
  alg.displayPaperStats(std::cout);
  {
    std::ofstream statsOut(mParams.inputFileNameStem(0) + ".stats");
    alg.displayStats(statsOut);
    alg.displayPaperStats(statsOut);
  }

  if (mGBParams.mOutputResult.value()) {
    // print basis
    {
      std::ofstream ogb(mParams.inputFileNameStem(0) + ".gb");
      ogb << "-- gb: ----\n";
      alg.getGB()->display(ogb);
    }
      
    // print syzygy basis
    {
      std::ofstream syzygyOut(mParams.inputFileNameStem(0) + ".syz");
      syzygyOut << "-- syz: ----\n";
      alg.getSyzTable()->display(syzygyOut);
      syzygyOut << std::endl;
    }
  }
}

const char* SigGBAction::staticName() {
  return "siggb";
}

const char* SigGBAction::name() const {
  return staticName();
}

const char* SigGBAction::description() const {
  return "Compute a signature Grobner basis. "
    "The project name is an optional direct parameter.";
}

const char* SigGBAction::shortDescription() const {
  return "Compute a signature Grobner basis";
}

void SigGBAction::pushBackParameters(
  std::vector<mic::CliParameter*>& parameters
) {
  mParams.pushBackParameters(parameters);
  mGBParams.pushBackParameters(parameters);
  parameters.push_back(&mUseSingularCriterionEarly);
  parameters.push_back(&mPostponeKoszul);
  parameters.push_back(&mUseBaseDivisors);
}

MATHICGB_NAMESPACE_END
