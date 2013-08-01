// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "mathicgb/stdinc.h"

#include "GBAction.hpp"
#include "SigGBAction.hpp"
#include "MatrixAction.hpp"
#include "HelpAction.hpp"
#include "mathicgb/LogDomainSet.hpp"
#include <mathic.h>
#include <cctype>
#include <iostream>
#include <exception>

// This is to satisfy the code checker that requires every file to contain
// these macroes.
MATHICGB_NAMESPACE_BEGIN
MATHICGB_NAMESPACE_END

int main(int argc, char **argv) {
  try {
    mathic::CliParser parser;
    parser.registerAction<mgb::SigGBAction>();
    parser.registerAction<mgb::GBAction>();
    parser.registerAction<mgb::MatrixAction>();
    parser.registerAction<mgb::HelpAction>();

    std::vector<std::string> commandLine(argv, argv + argc);
    commandLine.erase(commandLine.begin());

    parser.parse(commandLine)->performAction();
  } catch (const mathic::MathicException& e) {
    mathic::display(e.what());
    return -1;
  } catch (std::exception& e) {
    mathic::display(e.what());
    return -1;  
  } catch (...) {
    std::cout << "UNKNOWN ERROR" << std::endl;
    // maybe there is some outer exception handler that might say something
    // reasonable about this exception, so rethrow the exception.
    throw;
  }

  mgb::LogDomainSet::singleton().printReport(std::cerr);
  return 0;
};
