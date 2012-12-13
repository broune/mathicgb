#include "stdinc.h"
#include "LogDomain.hpp"

#include "LogDomainSet.hpp"
#include <iostream>

LogDomain<true>::LogDomain(
  const char* const name,
  const char* const description,
  const bool enabled
):
  mEnabled(enabled),
  mName(name),
  mDescription(description),
  mInterval()
{
  LogDomainSet::singleton().registerLogDomain(*this);
}

LogDomain<true>::~LogDomain() {
  if (enabled() && mHasTime) {
    stream() << mName << " total time:           ";
    mInterval.print(stream());
    stream() << '\n';
  }
}

std::ostream& LogDomain<true>::stream() {
  return std::cerr;
}

LogDomain<true>::Timer LogDomain<true>::timer() {
  return Timer(*this);
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

  MATHICGB_ASSERT(mName != 0);
  stream() << mName << " time recorded:        ";
  interval.print(stream());
  stream() << std::endl;
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
  mTimerRunning = false;
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

