#ifndef MATHICGB_NON_COPYABLE_GUARD
#define MATHICGB_NON_COPYABLE_GUARD

/// Derive from this class to disable the compiler-generated copy
/// constructor and assignment. T should be the class that is deriving
/// from NonCopyable.
///
/// The purpose of the template parameter is to avoid any chance of
/// getting a diamond-graph inheritance graph. Diamond graphs can lead
/// to runtime overhead.
template<class T>
class NonCopyable {
public:
  NonCopyable() {}
  NonCopyable(NonCopyable&&) {} // still movable.

private:
  NonCopyable(const NonCopyable&); // unavailable
  void operator=(const NonCopyable&); // unavailable
};

#endif
