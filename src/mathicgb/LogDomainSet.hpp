#ifndef MATHICGB_LOG_DOMAIN_SET_GUARD
#define MATHICGB_LOG_DOMAIN_SET_GUARD

#include "LogDomain.hpp"
#include <vector>
#include <algorithm>
#include <cstring>
#include <ostream>
#include <tbb/tbb.h>

class LogDomainSet {
public:
  void registerLogDomain(LogDomain<true>& domain);
  void registerLogDomain(const LogDomain<false>& domain) {}

  LogDomain<true>* logDomain(const char* const name);

  const std::vector<LogDomain<true>*>& logDomains() const {return mLogDomains;}

  void printReport(std::ostream& out) const;
  void printTimeReport(std::ostream& out) const;
  void printCountReport(std::ostream& out) const;

  static LogDomainSet& singleton();

private:
  LogDomainSet(); // private for singleton

  std::vector<LogDomain<true>*> mLogDomains;
  tbb::tick_count mStartTime;
};

#endif
