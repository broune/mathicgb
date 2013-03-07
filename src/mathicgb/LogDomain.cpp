#include "stdinc.h"
#include "LogDomain.hpp"

#include "LogDomainSet.hpp"
#include <mathic.h>
#include <tbb/tbb.h>
#include <iostream>

static const auto logDomainGlobalStartTime= tbb::tick_count::now();

LogDomain<true>::LogDomain(
  const char* const name,
  const char* const description,
  const bool enabled,
  const bool streamEnabled
):
  mEnabled(enabled),
  mStreamEnabled(streamEnabled),
  mName(name),
  mDescription(description),
  mInterval()
{
  LogDomainSet::singleton().registerLogDomain(*this);
}

std::ostream& LogDomain<true>::stream() {
  return std::cerr;
}

LogDomain<true>::Timer LogDomain<true>::timer() {
  return Timer(*this);
}

double LogDomain<true>::loggedSecondsReal() const {
  return mInterval.realSeconds;
}

void LogDomain<true>::TimeInterval::print(std::ostream& out) const {
  const auto oldFlags = out.flags();
  const auto oldPrecision = out.precision();
  out.precision(3);
  out << std::fixed << realSeconds << "s (real)";
  // todo: restore the stream state using RAII, since the above code might
  // throw an exception.
  out.precision(oldPrecision);
  out.flags(oldFlags);
}

void LogDomain<true>::recordTime(TimeInterval interval) {
  if (!enabled())
    return;
  mInterval.realSeconds += interval.realSeconds;
  mHasTime = true;

  if (streamEnabled()) {
    MATHICGB_ASSERT(mName != 0);
    stream() << mName << " time recorded:        ";
    interval.print(stream());
    stream() << std::endl;
  }
}

LogDomain<true>::Timer::Timer(LogDomain<true>& logger):
  mLogger(logger),
  mTimerRunning(false),
  mRealTicks()
{
  start();
}

LogDomain<true>::Timer::~Timer() {
  stop();
}

void LogDomain<true>::Timer::stop() {
  if (!running())
    return;
  mTimerRunning = false;
  if (!mLogger.enabled())
    return;
  TimeInterval interval;
  interval.realSeconds = (tbb::tick_count::now() - mRealTicks).seconds();
  mLogger.recordTime(interval);
  return;
}

void LogDomain<true>::Timer::start() {
  if (!mLogger.enabled() || mTimerRunning)
    return;
  mTimerRunning = true;
  mRealTicks = tbb::tick_count::now();
}

LogDomainInternal::LogAliasRegisterer::LogAliasRegisterer(const char* alias, const char* of) {
  LogDomainSet::singleton().registerLogAlias(alias, of);
}
