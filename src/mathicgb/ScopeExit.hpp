#ifndef MATHICGB_SCOPE_EXIT_GUARD
#define MATHICGB_SCOPE_EXIT_GUARD

// Guard holds an action to call and calls it unless it has been released.
template<class T>
class Guard {
public:
  ~Guard() {
    if (mOwning && mActive)
      mAction();
  }

private:
  friend struct GuardMaker;
  Guard(T&& action, const bool& active):
    mAction(std::move(action)), mOwning(true), mActive(active) {}

  // Most compilers should elide the call to this construtor, but it must be
  // here anyway and we should support even a crazy compiler that decides to
  // call it.
  Guard(Guard<T>&& guard):
    mAction(std::move(guard.mAction)), mOwning(true), mActive(guard.mActive)
  {
    assert(guard.mActive);
    guard.mOwning = false; // to avoid calling mAction twice
  }

  bool mOwning;
  const bool& mActive;
  const T mAction;
};

// The class user code interacts with to dismiss an action.
class Dismisser {
public:
  Dismisser(bool& active): mActive(active) {}
  void dismiss() {mActive = false;}

private:
  bool& mActive;
};

// Helper class that allows convenient syntax for the macro by overloading
// operator+.
struct GuardMaker {
public:
  GuardMaker(const bool& active): mActive(active) {}

  template<class T>
  Guard<T> operator+(T&& t) {return Guard<T>(std::forward<T>(t), mActive);}

private:
  const bool& mActive;
};

#define MYLIB__CAT_HELPER(A, B) A##B
#define MYLIB__CAT(A, B) MYLIB__CAT_HELPER(A, B)
#define MYLIB__UNIQUE(NAME) MYLIB__CAT(MyLib_,MYLIB__CAT(NAME,__LINE__))

// Example, with no need to dismiss:
//   FILE* file = fopen("file.txt", "r");
//   MATHICGB_SCOPE_EXIT() {
//     fclose(file);
//     std::cout << "file closed";
//   };
//   // ...
//   return; // the file is closed
//
// Example, with need to dismiss:
//   v.push_back(5);
//   MATHICGB_SCOPE_EXIT(name) {v.pop_back();};
//   // ...
//   if (error)
//     return; // the pop_back is done
//   name.dismiss();
//   return; // the pop_back is not done
//
// The middle line is a no-op if the name parameter expands to nothing.
// When NAME expands to nothing, we need the static_cast to prevent
// the compiler from parsing the middle line as a redeclaration of
// myBool. In the final line we use that a const reference keeps
// temporary objects alive until the end of the scope. It would be correct
// to copy, but this way we can keep the copy constructor private
// and help out any compiler that has issue with eliding that copy.
#define MATHICGB_SCOPE_EXIT(NAME) \
  bool MYLIB__UNIQUE(active) = true; \
  ::Dismisser NAME(static_cast<bool&>(MYLIB__UNIQUE(active))); \
  const auto& MYLIB__UNIQUE(guard) = \
    ::GuardMaker(MYLIB__UNIQUE(active)) + [&]

// Without this pragma, MSVC will say
//  warning C4003: not enough actual parameters for macro 'MYLIB_SCOPE_EXIT'
// when using MYLIB_SCOPE_EXIT without a name parameter. Not very happy about
// turning off the warning. I wonder if there is a way to avoid the warning in
// this case without turning it off.
#pragma warning (disable: 4003)

#endif
