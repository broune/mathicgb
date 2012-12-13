#ifndef MATHICGB_LOG_DOMAIN_GUARD
#define MATHICGB_LOG_DOMAIN_GUARD

#include <ostream>

template<bool Enabled>
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

  const char* name() const {return mName;}
  const char* description() const {return mDescription;}
  bool enabled() const {return mEnabled;}

  void setEnabled(const bool enabled) {mEnabled = enabled;}

  std::ostream& stream();

private:
  const char* mName;
  const char* mDescription;
  bool mEnabled;
};

template<>
class LogDomain<false> {
public:
  static const bool compileTimeEnabled = false;

  LogDomain(const char* const, const char* const, const bool) {}

  bool enabled() const {return false;}

  std::ostream& stream() {
    MATHICGB_ASSERT(false);
    return *static_cast<std::ostream*>(0);
  }
};

namespace LogDomainInternal {
  template<class Tag, bool Default>
  struct SelectValue {static const bool value = Default;};

  template<class> struct Tag_ {};
  template<class> struct Tag_0 {};
  template<class> struct Tag_1 {};

  template<bool Default>
  struct SelectValue<Tag_0<int>, Default> {static const bool value = false;};

  template<bool Default>
  struct SelectValue<Tag_1<int>, Default> {static const bool value = true;};
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

/// Defines logs::obj##NAME to be a LogDomain that is compile-time
/// enabled depending on MATHICGB_LOG_##NAME (see MATHICGB_CAPTURE_LOG_ENABLED)
/// and that is initially runtime enabled depending on the value of
/// DEFAULT_RUNTIME_ENABLED.
#define MATHICGB_DEFINE_LOG_DOMAIN_WITH_DEFAULTS(NAME, DESCRIPTION, DEFAULT_RUNTIME_ENABLED, DEFAULT_COMPILE_TIME_ENABLED) \
  MATHICGB_CAPTURE_LOG_ENABLED(NAME, DEFAULT_COMPILE_TIME_ENABLED); \
  namespace logs { \
    typedef LogDomain< ::LogDomainInternal::value_##NAME> Type##NAME; \
    Type##NAME NAME(#NAME, DESCRIPTION, DEFAULT_RUNTIME_ENABLED); \
  }

/// Defines logs::##NAME to be a LogDomain. It is compile-time enabled
/// by default and runtime disabled by default.
#define MATHICGB_DEFINE_LOG_DOMAIN(NAME, DESCRIPTION) \
  MATHICGB_DEFINE_LOG_DOMAIN_WITH_DEFAULTS(NAME, DESCRIPTION, 0, 1);

/// Have the code X in the program if and only if logging is not globally
/// compile-time disabled.
#define MATHICGB_IF_LOG(X) X

#define MATHICGB_LOG(DOMAIN) \
  if (::logs::##NAME.enabled()) ::logs::##NAME.stream()

#endif
