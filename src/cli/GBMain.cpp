#include "mathicgb/stdinc.h"

#include "GBAction.hpp"
#include "SigGBAction.hpp"
#include "MatrixAction.hpp"
#include <mathic.h>
#include <cctype>
#include <iostream>
#include <exception>

int main(int argc, char **argv) {
  try {
    mathic::CliParser parser;
    parser.registerAction<SigGBAction>();
    parser.registerAction<GBAction>();
    parser.registerAction<MatrixAction>();
    parser.registerAction<mathic::HelpAction>();

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
  return 0;
};

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
