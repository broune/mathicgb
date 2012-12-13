#include "stdinc.h"
#include "LogDomain.hpp"

#include "LogDomainSet.hpp"
#include <iostream>

LogDomain<true>::LogDomain(
  const char* const name,
  const char* const description,
  const bool enabled
):
  mName(name),
  mDescription(description),
  mEnabled(enabled)
{
  LogDomainSet::singleton().registerLogDomain(*this);
}

std::ostream& LogDomain<true>::stream() {
  return std::cerr;
}
