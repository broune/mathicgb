#include "mathicgb/stdinc.h"
#include "SigGBAction.hpp"

#include "mathicgb/Ideal.hpp"
#include "mathicgb/SignatureGB.hpp"
#include "mathicgb/io-util.hpp"
#include <fstream>
#include <iostream>

SigGBAction::SigGBAction():
  mUseSingularCriterionEarly("earlySingularCriterion",
    "Apply the singular S-pair elimination criterion before queueing "
    "that S-pair. Otherwise, the criterion is only checked just before "
    "the S-pair would otherwise cause a polynomial reduction to occur. "
    "This criterion is only relevant to the signature Buchberger "
    "algorithm.",
    false),

  mPostponeKoszul("postponeKoszul",
    "Postpone the construction of Koszul syzygy signatures.",
    true),

  mUseBaseDivisors("useBaseDivisors",
    "Use high ratio and low ratio base divisors to eliminate "
    "S-spairs quickly based on signature.",
    true),

  mModuleOrder("moduleOrder",
    "The free module term order.\n",
    4),

  mParams(1, 1)
{
  std::ostringstream orderOut;
  FreeModuleOrder::displayOrderTypes(orderOut);
  mModuleOrder.appendToDescription(orderOut.str());
}

void SigGBAction::directOptions(
  std::vector<std::string> tokens,
  mic::CliParser& parser
) {
  mParams.directOptions(tokens, parser);
}

void SigGBAction::performAction() {
  mParams.perform();
  mGBParams.perform();

  // read input file
  std::unique_ptr<Ideal> ideal;
  {
    const std::string inputIdealFile = mParams.inputFileNameStem(0) + ".ideal";
    std::ifstream inputFile(inputIdealFile.c_str());
    if (inputFile.fail())
      mic::reportError("Could not read input file \"" + inputIdealFile + '\n');
    ideal = Ideal::parse(inputFile);
  }
  std::unique_ptr<PolyRing const> ring(&(ideal->ring()));

  SignatureGB alg(
    std::move(*ideal),
    mModuleOrder.value(),
    Reducer::reducerType(mGBParams.mReducer.value()),
    mGBParams.mDivisorLookup.value(),
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
    std::ofstream statsOut((mParams.inputFileNameStem(0) + ".stats").c_str());
    alg.displayStats(statsOut);
    alg.displayPaperStats(statsOut);
  }

  if (mGBParams.mOutputResult.value())
    {
      // print basis
      {
        std::ofstream ogb((mParams.inputFileNameStem(0) + ".gb").c_str());
        ogb << "-- gb: ----\n";
        alg.getGB()->display(ogb);
      }
      
      // print syzygy basis
      {
        std::ofstream syzygyOut((mParams.inputFileNameStem(0) + ".syz").c_str());
        syzygyOut << "-- syz: ----\n";
        alg.getSyzTable()->display(syzygyOut, 1);
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
  parameters.push_back(&mModuleOrder);
}
