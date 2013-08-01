// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "mathicgb/stdinc.h"
#include "HelpAction.hpp"

#include "mathicgb/LogDomain.hpp"
#include "mathicgb/LogDomainSet.hpp"
#include <iostream>

MATHICGB_NAMESPACE_BEGIN

void HelpAction::performAction() {
  if (topic() != "logs") {
    mathic::HelpAction::performAction();
    return;
  }

  const char* header =
    "MathicGB offers a set of logs that can be enabled or disabled "
    "individually.\n"
    "\n"
    "A log can have a streaming component and a summary "
    "component. The streaming component is printed as the logged event "
    "occurs while the summary component is printed at the end. For example "
    "a log for opening files would stream a message every time a file is "
    "opened and it would report the total number of files opened when the "
    "program has run to the end. The summary can also include how long "
    "something took.\n"
    "\n"
    "An enabled log always prints its summary if it registered any events. "
    "The streaming component can be turned off even for an enabled log, and "
    "this is the default for some logs.\n"
    "\n"
    "You specify the log configuration using the option -log X, where X is "
    "a comma-seperated list of log specifications. Here is an example:\n"
    "\n"
    "    A,+B,-C,D+,E-\n"
    "\n"
    "This enables logs A, B, D and E and disables C. Furthermore, streaming "
    "for D is turned on while it is turned off for E. The streaming setting "
    "for A, B and C is the default for those logs.\n"
    "\n"
    "A prefix of - disables the log while no prefix or a prefix of + enables "
    "it. A suffix of - turns off streaming while a suffix of + turns it on. "
    "If there is no suffix then the setting for streaming is unchanged. A "
    "prefix or suffix of 0 means do nothing.\n"
    "\n"
    "The following is a list of all compile-time enabled logs. The prefixes "
    "and suffixes indicate the default state of the log.\n";
  mathic::display(header);
  auto& logs = LogDomainSet::singleton().logDomains();
  for (auto it = logs.begin(); it != logs.end(); ++it) {
    const auto toSign = [](const bool b) {return b ? '+' : '-';};
    std::cerr
      << "\n "
      << toSign((*it)->enabled())
      << (*it)->name()
      << toSign((*it)->streamEnabledPure())
      << '\n';
    mathic::display((*it)->description(), "   ");
  }

  const char* aliasDescription =
    "\nA log alias is a short-hand where one name stands for several log "
    "specifications. Prefixes and Suffixes also apply to aliases. If X "
    "expands to A+,-B,C then +X- expands to +A-,+B-,+C-. 0all- turns off "
    "all streaming output without enabling or disabling any logs.\n"
    "\n"
    "The following is a list of all aliases.\n";
  mathic::display(aliasDescription);
  auto& aliases = LogDomainSet::singleton().aliases();
  for (auto it = aliases.begin(); it != aliases.end(); ++it) {
    std::cerr << "\n " << it->first << " expands to\n";
    std::string str = it->second;
    std::replace(str.begin(), str.end(), ',', ' ');
    mathic::display(str, "   ");
  }
  std::cerr <<
    "\n none expands to nothing\n"
    "\n"
    " all expands to all log names\n";
    "\n";
}

MATHICGB_NAMESPACE_END
