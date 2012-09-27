#include "mathicgb/stdinc.h"

#include "mathicgb/Ideal.hpp"
#include "mathicgb/SignatureGB.hpp"
#include "mathicgb/BuchbergerAlg.hpp"
#include "mathicgb/io-util.hpp"
#include "mathicgb/MTArray.hpp"

#include <mathic.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <istream>
#include <cctype>

extern int tracingLevel;

class CliActionSignature : public mic::Action {
public:
  CliActionSignature():
    mUseSingularCriterionEarly("earlySingularCriterion",
      "Apply the singular S-pair elimination criterion before queueing "
      "that S-pair. Otherwise, the criterion is only checked just before "
      "the S-pair would otherwise cause a polynomial reduction to occur. "
      "This criterion is only relevant to the signature Buchberger "
      "algorithm.",
      false),

    mPreferSparseReducers("preferSparseReducers",
      "If true, always use the sparsest reducer in polynomial reduction. "
      "This option impacts both classic and signature constrained "
      "polynomial reduction. Ties are broken by taking the oldest reducer. "
      "If this option is false, the oldest reducer is always used.",
      true),

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

    mClassicBuchbergerAlgorithm("classicBuchbergerAlgorithm",
      "Use the classic buchberger algorithm. If true, ignore "
      "any options having to do with signatures or the signature "
      "algorithm.",
      false),

    mPostponeKoszul("postponeKoszul",
      "Postpone the construction of Koszul syzygy signatures.",
      true),

    mUseBaseDivisors("useBaseDivisors",
      "Use high ratio and low ratio base divisors to eliminate "
      "S-spairs quickly based on signature.",
      true),

    mSPairQueue("spairQueue",
      "The priority queue used to order S-pairs.\n"
      "  0   tournament tree in front of triangle\n"
      "  1   heap in front of triangle\n"
      "  2   tournament tree\n"
      "  3   heap\n",
      0),

    mTracingLevel("tracingLevel",
      "How much information to print out about actions taken. Nothing "
      "extra is printed if the value is zero. Higher values "
      "result in more information.",
      0),

    mBreakAfter("breakAfter",
      "Stop the computation after this many elements have been added to "
      "the basis. The computation runs uninterrupted if the value is zero.",
      0),

    mPrintInterval("printInterval",
      "Print information about the computation every time this many S-pair "
      "reductions have been performed. Do not print information like this "
      "during the computation if the value is zero.",
      0),

    mMonomialTable("monomialTable",
      "The kind of monomial table data structure to use.\n",
      2),

    mDivisorLookup("divisorLookup",
      "The divisor lookup data structure to use.\n",
      2),

    mReducer("reducer",
      "The data structure to use for polynomial reduction.\n",
      4),

    mModuleOrder("moduleOrder",
      "The free module term order.\n",
      4),

    mProjectName("projectName",
      "For a value X, read the input ideal from X.ideal, "
      "write statistics to X.stats and write the output to X.gb",
      ""),

    mStrategy("strategy",
      "The integer representing the strategy to use for the computation:\n"
      "  0   signature+gbstop\n"
      "  1   buchberger\n"
      "  2   late koszul\n"
      "  4   not used (used to be MES reduction)\n"
      "  8   signature\n",
      2)
  {
    {
      std::ostringstream orderOut;
      FreeModuleOrder::displayOrderTypes(orderOut);
      mModuleOrder.appendToDescription(orderOut.str());
    }
    {
      std::ostringstream reducerOut;
      Reducer::displayReducerTypes(reducerOut);
      mReducer.appendToDescription(reducerOut.str());
    }
    {
      std::ostringstream divisorLookupOut;
      DivisorLookup::displayDivisorLookupTypes(divisorLookupOut);
      mDivisorLookup.appendToDescription(divisorLookupOut.str());
    }
    {
      std::ostringstream monomialTableOut;
      MonomialTableArray::displayMTTypes(monomialTableOut);
      mMonomialTable.appendToDescription(monomialTableOut.str());
    }
  }

  virtual ~CliActionSignature() {}

  virtual void directOptions(
    std::vector<std::string> tokens,
    mic::CliParser& parser
  ) {
    switch (tokens.size()) {
    case 8: mStrategy.processArgument(tokens[6]);
      mClassicBuchbergerAlgorithm.setValue((mStrategy.value() &  1) == 1);
      mPostponeKoszul.setValue((mStrategy.value() &  2) != 0);
      mUseBaseDivisors.setValue((mStrategy.value() & 16) != 0);
    case 7: mTracingLevel.processArgument(tokens[5]);
    case 6: // ignore this value
    case 5: mMonomialTable.processArgument(tokens[3]);
    case 4: mDivisorLookup.processArgument(tokens[2]);
    case 3: mReducer.processArgument(tokens[1]);
    case 2: mModuleOrder.processArgument(tokens[0]);
    case 1: mProjectName.processArgument(tokens.back());
    case 0: break;
    default: mic::reportError
      (std::string("Too many direct options for action ") + staticName());
    }
  }

  virtual void performAction() {
    tracingLevel = mTracingLevel.value();

    // read input file
    std::unique_ptr<Ideal> ideal;
    {
      std::string const inputIdealFile = mProjectName.value() + ".ideal";
      std::ifstream inputFile(inputIdealFile.c_str());
      if (inputFile.fail())
        mic::reportError("Could not read input file " + inputIdealFile);
      ideal.reset(Ideal::parse(inputFile));
    }
    std::unique_ptr<PolyRing const> ring(&(ideal->ring()));
      
    if (mClassicBuchbergerAlgorithm.value()) {
      BuchbergerAlg alg(
        *ideal,
        mModuleOrder.value(),
        Reducer::ReducerType(mReducer.value()),
        mDivisorLookup.value(),
        mPreferSparseReducers.value(),
        mSPairQueue.value());
      alg.setBreakAfter(mBreakAfter.value());
      alg.setPrintInterval(mPrintInterval.value());
      alg.setUseAutoTopReduction(mAutoTopReduce.value());
      alg.setUseAutoTailReduction(mAutoTailReduce.value());

      alg.computeGrobnerBasis();
      alg.printStats(std::cerr);

      std::ofstream statsOut((mProjectName.value() + ".stats").c_str());
      alg.printStats(statsOut);
      std::ofstream gbOut((mProjectName.value() + ".gb").c_str());
      output(gbOut, alg.basis());
    } else {
      SignatureGB alg(
        *ideal,
        mModuleOrder.value(),
        Reducer::reducerType(mReducer.value()),
        mDivisorLookup.value(),
        mMonomialTable.value(),
        mPostponeKoszul.value(),
        mUseBaseDivisors.value(),
        mPreferSparseReducers.value(),
        mUseSingularCriterionEarly.value(),
        mSPairQueue.value());
      alg.setBreakAfter(mBreakAfter.value());
      alg.setPrintInterval(mPrintInterval.value());
      alg.computeGrobnerBasis();

      // print statistics
      alg.displayStats(std::cout);
      alg.displayPaperStats(std::cout);

      {
        std::ofstream statsOut((mProjectName.value() + ".stats").c_str());
        alg.displayStats(statsOut);
        alg.displayPaperStats(statsOut);
      }

      // print basis
      {
        std::ofstream ogb((mProjectName.value() + ".gb").c_str());
        ogb << "-- gb: ----\n";
        alg.getGB()->display(ogb);
      }

      // print syzygy basis
      {
        std::ofstream syzygyOut((mProjectName.value() + ".syz").c_str());
        syzygyOut << "-- syz: ----\n";
        alg.getSyzTable()->display(syzygyOut, 1);
        syzygyOut << std::endl;
      }
    }
  }

  static const char* staticName() {return "gb";}

  virtual const char* name() const {return staticName();}
  virtual const char* description() const {
    return
      "Compute a Grobner basis. Takes optional direct parameters in this order:\n"
      "  <module order> <reducer> <divisor lookup> <monomial table>\n"
      "  <ignore this value> <trace level> <strategy> <project name>\n"
      "If fewer than 8 direct parameters, only so many options are read. They "
      "are read in the order given, except that the the last direct "
      "parameter, if any, is always assumed to be the project name. "
      "Parameters given with a dash take priority over direct parameters.";
  }
  virtual const char* shortDescription() const {return "Compute a Grobner basis";}
  
  virtual void pushBackParameters(std::vector<mic::CliParameter*>& parameters) {
    parameters.push_back(&mUseSingularCriterionEarly);
    parameters.push_back(&mPreferSparseReducers);
    parameters.push_back(&mAutoTailReduce);
    parameters.push_back(&mAutoTopReduce);
    parameters.push_back(&mClassicBuchbergerAlgorithm);
    parameters.push_back(&mPostponeKoszul);
    parameters.push_back(&mUseBaseDivisors);
    parameters.push_back(&mSPairQueue);
    parameters.push_back(&mTracingLevel);
    parameters.push_back(&mBreakAfter);
    parameters.push_back(&mPrintInterval);
    parameters.push_back(&mMonomialTable);
    parameters.push_back(&mDivisorLookup);
    parameters.push_back(&mReducer);
    parameters.push_back(&mModuleOrder);
    parameters.push_back(&mProjectName);

    // do not expose the strategy parameter - it is only here to support
    // the old format of using direct numeric parameters to the action.
    //parameters.push_back(&mStrategy);
  }

private:
  mic::BoolParameter mUseSingularCriterionEarly;
  mic::BoolParameter mPreferSparseReducers;
  mic::BoolParameter mAutoTailReduce;
  mic::BoolParameter mAutoTopReduce;
  mic::BoolParameter mClassicBuchbergerAlgorithm;
  mic::BoolParameter mPostponeKoszul;
  mic::BoolParameter mUseBaseDivisors;
  mic::IntegerParameter mSPairQueue;
  mic::IntegerParameter mTracingLevel;
  mic::IntegerParameter mBreakAfter;
  mic::IntegerParameter mPrintInterval;
  mic::IntegerParameter mMonomialTable;
  mic::IntegerParameter mDivisorLookup;
  mic::IntegerParameter mReducer;
  mic::IntegerParameter mModuleOrder;
  mic::StringParameter mProjectName;
  mic::IntegerParameter mStrategy;
};

int oldmain(int argc, char **argv);
int main(int argc, char **argv) {
  //oldmain(argc, argv);
  try {
    mic::CliParser parser;
    parser.registerAction<CliActionSignature>();
    parser.registerAction<mic::HelpAction>();
    std::vector<std::string> commandLine(argv, argv + argc);
    if (commandLine.size() >= 2 &&
      !commandLine[1].empty() &&
      std::isdigit(commandLine[1][0])) {
      commandLine[0] = "gb";
    } else
      commandLine.erase(commandLine.begin());
    // todo: remove the .release once parser returns unique_ptr
    // instead of auto_ptr.
    std::unique_ptr<mic::Action> action(parser.parse(commandLine).release());
    action->performAction();
  } catch (const mic::MathicException& e) {
    mic::display(e.what());
    return -1;
  }
  return 0;
};

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
