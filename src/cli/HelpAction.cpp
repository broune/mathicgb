#include "mathicgb/stdinc.h"
#include "HelpAction.hpp"

#include "mathicgb/LogDomain.hpp"
#include "mathicgb/LogDomainSet.hpp"

#include <iostream>

void HelpAction::performAction() {
  if (topic() != "logs") {
    mathic::HelpAction::performAction();
    return;
  }

  auto& logs = LogDomainSet::singleton().logDomains();
  for (auto it = logs.begin(); it != logs.end(); ++it) {
    std::cout << "\n " << (*it)->name() << '\n';
    mathic::display((*it)->description(), "   ");
  }
}
