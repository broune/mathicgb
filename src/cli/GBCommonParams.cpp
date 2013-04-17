#include "mathicgb/stdinc.h"
#include "GBCommonParams.hpp"

#include "mathicgb/MTArray.hpp"
#include "mathicgb/PolyReducer.hpp"
#include "mathicgb/DivisorLookup.hpp"

GBCommonParams::GBCommonParams():
  mPreferSparseReducers("preferSparseReducers",
    "If true, always use the sparsest reducer in polynomial reduction. "
    "This option impacts both classic and signature constrained "
    "polynomial reduction. Ties are broken by taking the oldest reducer. "
    "If this option is false, the oldest reducer is always used.",
    true),

  mOutputResult("outputResult",
    "If true, output the resulting Groebner or signature basis "
    "to the file <projectName>.gb and in the signature basis "
    "case, the signatures of the syzygies are placed in <projectName>.syz",
    false),

  mSPairQueue("spairQueue",
    "The priority queue used to order S-pairs.\n"
    "  0   tournament tree in front of triangle\n"
    "  1   heap in front of triangle\n"
    "  2   tournament tree\n"
    "  3   heap\n",
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

  mMemoryQuantum("memoryQuantumForReducer",
    "Specifies how many items to allocate memory for at a time for the reducer.",
    1024 * 1024)
{
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

void GBCommonParams::pushBackParameters(
  std::vector<mathic::CliParameter*>& parameters
) {
  parameters.push_back(&mPreferSparseReducers);
  parameters.push_back(&mOutputResult);
  parameters.push_back(&mSPairQueue);
  parameters.push_back(&mBreakAfter);
  parameters.push_back(&mPrintInterval);
  parameters.push_back(&mMonomialTable);
  parameters.push_back(&mDivisorLookup);
  parameters.push_back(&mReducer);
  parameters.push_back(&mMemoryQuantum);
}

void GBCommonParams::perform() {
  // currently there is nothing to do
}
