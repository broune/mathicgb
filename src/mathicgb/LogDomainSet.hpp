#ifndef MATHICGB_LOG_DOMAIN_SET_GUARD
#define MATHICGB_LOG_DOMAIN_SET_GUARD

#include "LogDomain.hpp"
#include <vector>
#include <algorithm>

class LogDomainSet {
public:
  void registerLogDomain(LogDomain<true>& domain) {
    mLogDomains.push_back(&domain);
  }

  void registerLogDomain(const LogDomain<false>& domain) {
  }

  LogDomain<true>* logDomain(const char* const name) {
    const auto func = [&](const LogDomain<true>* const ld){
      return std::strcmp(ld->name(), name) == 0;
    };
    const auto it = std::find_if(mLogDomains.begin(), mLogDomains.end(), func);
    return it == mLogDomains.end() ? static_cast<LogDomain<true>*>(0) : *it;
  }

  const std::vector<LogDomain<true>*> logDomains() const {return mLogDomains;}

  static LogDomainSet& singleton() {
    static LogDomainSet set;
    return set;
  }

private:
  LogDomainSet() {}

  std::vector<LogDomain<true>*> mLogDomains;
};


#endif
