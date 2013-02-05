#include "stdinc.h"
#include "LogDomainSet.hpp"

#include <mathic.h>

LogDomainSet::LogDomainSet():
  mStartTime(tbb::tick_count::now()) {
}

void LogDomainSet::registerLogDomain(LogDomain<true>& domain) {
  mLogDomains.push_back(&domain);
}

LogDomain<true>* LogDomainSet::logDomain(const char* const name) {
  const auto func = [&](const LogDomain<true>* const ld){
    return std::strcmp(ld->name(), name) == 0;
  };
  const auto it = std::find_if(mLogDomains.begin(), mLogDomains.end(), func);
  return it == mLogDomains.end() ? static_cast<LogDomain<true>*>(0) : *it;
}

void LogDomainSet::printReport(std::ostream& out) const {
  const auto allTime = (tbb::tick_count::now() - mStartTime).seconds();

  mathic::ColumnPrinter pr;
  auto& names = pr.addColumn(true);
  auto& times = pr.addColumn(false);
  auto& ratios = pr.addColumn(false);
  times.precision(3);
  times << std::fixed;
  ratios.precision(3);
  ratios << std::fixed;

  names << "Log name  \n";
  times << "  Time/s (real)\n";
  ratios << "  Ratio\n";
  pr.repeatToEndOfLine('-');

  double timeSum = 0;
  bool somethingToReport = false;
  const auto end = logDomains().cend();
  for (auto it = logDomains().cbegin(); it != end; ++it) {
    const auto& log = **it;
    if (!log.enabled() || !log.hasTime())
      continue;
    somethingToReport = true;

    const auto logTime = log.loggedSecondsReal();
    timeSum += logTime;
    names << log.name() << '\n';
    times << logTime << '\n';
    ratios << mathic::ColumnPrinter::percentDouble(logTime, allTime) << '\n';
  }
  if (!somethingToReport)
    return;
  pr.repeatToEndOfLine('-');
  names << "sum\n";
  times << timeSum;
  ratios << mathic::ColumnPrinter::percentDouble(timeSum, allTime) << '\n';

  const auto oldFlags = out.flags();
  const auto oldPrecision = out.precision();
  out << std::fixed;
  out.precision(3);
  out << "***** Logging report *****\nTime elapsed: "
    << allTime << "s\n\n" << pr << '\n';

  // todo: restore the stream state using RAII, since the above code might
  // throw an exception.
  out.precision(oldPrecision);
  out.flags(oldFlags);
}

LogDomainSet& LogDomainSet::singleton() {
  static LogDomainSet set;
  return set;
}
