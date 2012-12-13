#ifndef MATHICGB_LOG_DOMAIN_GUARD
#define MATHICGB_LOG_DOMAIN_GUARD

#include <tbb/tbb.h>
#include <ostream>
#include <ctime>
#include <sstream>

/// A named area of logging that can be turned on or off at runtime and at
/// compile time.
///
/// A logger that is turned off at compile time emits no code
/// into the executable and all the code that writes to that logger is also
/// removed by the optimizer if it is written in the correct way. Use the
/// logging macroes to ensure proper use so that compile-time disabled
/// LogDomains properly have zero overhead. LogDomains can be turned on
/// and off at compile time and at runtime individually.
///
/// Compile-time enabled loggers automatically register themselves with
/// LogDomainSet::singleton().
///
/// @todo: support turning all loggers off globally with a macro, regardless
/// of their individual compile-time on/off setting.

template<bool CompileTimeEnabled>
class LogDomain {};

template<>
class LogDomain<true> {
public:
  static const bool compileTimeEnabled = true;

  LogDomain(
    const char* const name,
    const char* const description,
    const bool enabled
  );
  ~LogDomain();

  const char* name() const {return mName;}
  const char* description() const {return mDescription;}
  bool enabled() const {return mEnabled;}

  void setEnabled(const bool enabled) {mEnabled = enabled;}

  std::ostream& stream();

  /// Class for recording time that is logged.
  class Timer;

  /// Returns a started timer that you can move from.
  Timer timer();

private:
  struct TimeInterval {
    // todo: support user time too. clock() doesn't seem to sum the time
    // for all threads, so that didn't work.
    double realSeconds;

    void print(std::ostream& out) const;
  };
  void recordTime(TimeInterval interval);

  bool mEnabled;
  const char* mName;
  const char* mDescription;

  TimeInterval mInterval; /// Total amount of time recorded on this log.

  /// Indicates if any period of time has been recorded, even if that period
  /// of time was recorded as 0 seconds.
  bool mHasTime;
};

class LogDomain<true>::Timer {
public:
  /// Start the timer running. The elapsed time will be logged to the logger
  /// once the timer is stopped or destructed.
  Timer(LogDomain<true>& logger);

  /// Stops the timer.
  ~Timer();

  /// Returns true if the timer is currently recording time.
  bool running() const {return mTimerRunning;}

  /// Stops recording time and logs the elapsed time to the logger.
  ///
  /// This is a no-op if the timer is not running. If the logger
  /// is disabled then no time is logged.
  void stop();

  /// Start recording time on a stopped timer.
  ///
  /// This is a no-op is the timer is already running or if the logger is
  /// disabled.
  void start();

private:
  LogDomain<true>& mLogger;
  bool mTimerRunning;
  tbb::tick_count mRealTicks; // high precision
};

/// This is a compile-time disabled logger.
template<>
class LogDomain<false> {
public:
  static const bool compileTimeEnabled = false;

  LogDomain(const char* const, const char* const, const bool) {}

  bool enabled() const {return false;}

  class Timer {
  public:
    Timer(LogDomain<false>&) {}
    bool running() const {return false;}
    void stop() {}
    void start() {}
  };
  Timer timer() {return Timer(*this);}

  std::ostream& stream() {
    MATHICGB_ASSERT(false);
    return *static_cast<std::ostream*>(0);
  }
};

namespace LogDomainInternal {
  // Support code for the logging macroes


  template<class Tag, bool Default>
  struct SelectValue {static const bool value = Default;};

  template<class> struct Tag_ {};
  template<class> struct Tag_0 {};
  template<class> struct Tag_1 {};

  template<bool Default>
  struct SelectValue<Tag_0<int>, Default> {static const bool value = false;};

  template<bool Default>
  struct SelectValue<Tag_1<int>, Default> {static const bool value = true;};

  template<class L>
  struct LambdaRunner {L& log;};
  template<class L>
  LambdaRunner<L> lambdaRunner(L& log) {return LambdaRunner{log};}
  template<class L, class T>
  void operator+(LambdaRunner<L> runner, T& lambda) {lambda(runner.log);}
}

/// Defines LogDomainInternal::value_##NAME to be equal to the value of
/// the macro MATHICGB_LOG_##NAME if that macro expands to 0 or 1. Otherwise
/// the macro MATHICGB_LOG_##NAME is ignored and instead DEFAULT_VALUE is used.
#define MATHICGB_CAPTURE_LOG_ENABLED(NAME, DEFAULT_VALUE) \
  namespace LogDomainInternal { \
    template<class> struct Tag_MATHICGB_LOG_##NAME {}; \
    typedef MATHICGB_CONCATENATE_AFTER_EXPANSION(Tag_, MATHICGB_LOG_##NAME)<int> \
      SelectedTag_##NAME; \
    static const bool value_##NAME = \
      SelectValue<SelectedTag_##NAME, DEFAULT_VALUE>::value; \
  }

/// Defines a LogDomain with the given name and description.
///
/// The logger is default compile-time enabled depending on MATHICGB_LOG_##NAME
/// (see MATHICGB_CAPTURE_LOG_ENABLED) and it is initially runtime
/// enabled depending on the value of DEFAULT_RUNTIME_ENABLED.
#define MATHICGB_DEFINE_LOG_DOMAIN_WITH_DEFAULTS(NAME, DESCRIPTION, DEFAULT_RUNTIME_ENABLED, DEFAULT_COMPILE_TIME_ENABLED) \
  MATHICGB_CAPTURE_LOG_ENABLED(NAME, DEFAULT_COMPILE_TIME_ENABLED); \
  namespace logs { \
    typedef LogDomain<::LogDomainInternal::value_##NAME> Type##NAME; \
    Type##NAME NAME(#NAME, DESCRIPTION, DEFAULT_RUNTIME_ENABLED); \
  }

/// Defines a LogDomain with the given name and description.
///
/// By default, the logger is compile-time enabled and runtime disabled.
#define MATHICGB_DEFINE_LOG_DOMAIN(NAME, DESCRIPTION) \
  MATHICGB_DEFINE_LOG_DOMAIN_WITH_DEFAULTS(NAME, DESCRIPTION, 0, 1);

/// This expression yields an l-value reference to the indicated logger.
///
/// Example:
///   auto timer = MATHICGB_LOGGER(MyDomain).timer();
#define MATHICGB_LOGGER(DOMAIN) ::logs::##DOMAIN

/// This expression yields the type of the indicated logger.
///
/// Example:
///   if (MATHICGB_LOGGER_TYPE(MyDomain)::compileTimeEnabled)
///     std::ostream << "MyDomain is compiled time enabled";
#define MATHICGB_LOGGER_TYPE(DOMAIN) ::logs::Type##DOMAIN

/// Runs the code in the following scope delimited by braces {} if the indicated
/// logger is enabled - otherwise does nothing. Within the following scope
/// there is a local reference variable log that refers to the indicated
/// logger.
///
/// Example:
///   MATHICGB_IF_LOG(MyDomain) {
///     std::string msg;
///     expensiveFunction(msg);
///     log << msg;
///   }
#define MATHICGB_IF_LOG(DOMAIN) \
  if (MATHICGB_LOGGER(DOMAIN).enabled()) \
    LogDomainInternal::lambdaRunner(MATHICGB_LOGGER(DOMAIN)) + \
      [&](MATHICGB_LOGGER_TYPE(DOMAIN)& log)

/// Display information to the log using <<.
/// If domain is not enabled then the log message is not displayed and the
/// code after << is not executed.
///
/// Example: (f() only called if logger is enabled)
///   MATHICGB_LOG(domain) << "f() = " << f();
#define MATHICGB_LOG(DOMAIN) \
  if (MATHICGB_LOGGER(DOMAIN).enabled()) MATHICGB_LOGGER(DOMAIN).stream()

/// Will log the time to execute the remaining code in the current scope
/// to the indicated domain. Also supports printing a message using <<.
/// The message is printed right away while the time is printed when
/// the scope ends.
///
/// Example:
///   MATHICGB_LOG_SCOPE_TIME(MyDomain) << "Starting timed task";
#define MATHICGB_LOG_TIME(DOMAIN) \
  auto MATHICGB_timer##DOMAIN##_##__LINE__##(MATHICGB_LOGGER(DOMAIN).timer()); \
  MATHICGB_LOG(DOMAIN)

#endif
