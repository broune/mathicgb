#ifndef MATHICGB_LOG_DOMAIN_SET_GUARD
#define MATHICGB_LOG_DOMAIN_SET_GUARD

#include "LogDomain.hpp"
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <ostream>
#include <tbb/tbb.h>

class LogDomainSet {
public:
  void registerLogDomain(LogDomain<true>& domain);
  void registerLogDomain(const LogDomain<false>& domain) {}

  /// A log command has the format AXB, where
  ///   X       the name of a compile-time enabled log domain
  ///   A       a prefix
  ///   B       a suffix
  /// The possible values of A are
  ///           enabled X (this is the empty string)
  ///   +       enabled X
  ///   -       disable X
  /// The possible values of B are
  ///           leaving sub-state as-is (this is the empty string)
  ///   +       stream-enabled X
  ///   -       stream-disable X
  ///
  /// No white-space is allowed.
  /// If the command cannot be parsed then you will get an exception.
  ///
  /// *** Example ***
  /// Consider this sequence of commands:
  ///   "+MyLog-" will enabled MyLog, but silence any streaming from it.
  ///   "+MyLog" will make no difference, as MyLog is already enabled.
  ///   "-MyLog+" will disable MyLog, but set the streaming state to enabled.
  ///     As MyLog is disabled there will still be no streaming output.
  ///   "+MyLog" will enabled MyLog. Since the streaming state was enabled
  ///     before, we now get streaming.
  ///
  void performLogCommand(std::string cmd);

  /// Performs a comma-seperated list of commands. No white-space is allowed.
  void performLogCommands(const std::string& cmds);

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
