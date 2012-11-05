#ifdef MATHICGB_STDINC_GUARD
#error stdinc.h included twice. Only include stdinc.h once per cpp file.
#endif
#define MATHICGB_STDINC_GUARD

#ifdef _MSC_VER // For Microsoft Compiler in Visual Studio C++.

// Sometimes you know that a function will be called very rarely so you want to
// tell the compiler not to inline it even if it could be inlined at only a
// modest increase in code size. That is what MATHICGB_NO_INLINE does.
#define MATHICGB_NO_INLINE __declspec(noinline)

// Sometimes the compiler just will not inline functions that should
// be inlined. Use sparingly --- preferably only if a profiler says
// that a tiny often-called function consumes a significant amount of time.
#define MATHICGB_INLINE __forceinline

// Tells the compiler to always assume that the expression X is true.
#define MATHICGB_ASSUME(X) __assume(X)

// As MATHICGB_ASSUME, but might actually evaluate X at run-time if it has
// side-effects. The point is that this can be used on compilers with no other
// support for assuming things. So there is no difference on MS VC++.
#define MATHICGB_ASSUME_AND_MAY_EVALUATE(X) __assume(X)

// Tells the compiler that this function returns a pointer that is not an alias
// for any other point that is currently valid in the program - like malloc.
#define MATHICGB_RETURN_NO_ALIAS __declspec(restrict)

// Tells the compiler that this function will never throw an exception.
#define MATHICGB_NOTHROW __declspec(nothrow)

// Tells the compiler that this function has no effects except the return value
// and the return value depends only on the arguments and first-level
// indirections of the arguments. (this is the common denominator of GCC
// and MS VC++ capabilities)
#define MATHICGB_PURE __declspec(noalias)

// Tells the compiler that the return value of this function must be looked
// at by the caller. For example this is appropriate for realloc.
#define MATHICGB_MUST_CHECK_RETURN_VALUE

// Tells the compiler that the current line of code cannot be reached.
#define MATHICGB_UNREACHABLE __assume(false)

// Tells the compiler that a variable that is a pointer (not a reference)
// does not alias any other pointer that is used in the current scope.
#define MATHICGB_RESTRICT __restrict


#pragma warning (disable: 4996) // std::copy on pointers is flagged as dangerous
#pragma warning (disable: 4290) // VC++ ignores throw () specification.
#pragma warning (disable: 4127) // Warns about using "while (true)".
#pragma warning (disable: 4100) // Warns about unused parameters.
#pragma warning (disable: 4800) // Warns on int to bool conversion.
#pragma warning (disable: 4146) // Warns on unary minus on unsigned (bit trick)

// This warning warns about using the this pointer in base member
// initializer lists. This is a pretty good warning as that can
// obviously easily go wrong, but it is pretty useful to do as well,
// so the warning is turned off.
#pragma warning (disable: 4355)

#if defined (_M_IX86) || defined(_M_X64) // if on x86 (32 bit) or x64 (64 bit)
#define MATHICGB_USE_CUSTOM_ATOMIC_X86_X64_MSVC_4BYTE
#endif
#ifdef _M_X64 // if on x64 (64 bit)
#define MATHICGB_USE_CUSTOM_ATOMIC_X86_X64_MSVC_8BYTE
#endif

#elif defined (__GNUC__) // GCC compiler

#define MATHICGB_NO_INLINE __attribute__((noinline))
#define MATHICGB_INLINE __attribute__((always_inline)) inline
#define MATHICGB_ASSUME(X)
#define MATHICGB_ASSUME_AND_MAY_EVALUATE(X) do {if(!(X)){MATHICGB_UNREACHABLE;}while(0)}
#define MATHICGB_RETURN_NO_ALIAS __attribute__(malloc)
#define MATHICGB_NOTHROW __attribute__(nothrow)
#define MATHICGB_PURE __attribute__(pure)
#define MATHICGB_MUST_CHECK_RETURN_VALUE __attribute__(warn_unused_result)
#define MATHICGB_UNREACHABLE __builtin_unreachable()
#define MATHICGB_RESTRICT __restrict

#else

#define MATHICGB_NO_INLINE
#define MATHICGB_INLINE inline
#define MATHICGB_ASSUME(X)
#define MATHICGB_ASSUME_AND_MAY_EVALUATE(X)
#define MATHICGB_RETURN_NO_ALIAS
#define MATHICGB_NOTHROW
#define MATHICGB_PURE
#define MATHICGB_MUST_CHECK_RETURN_VALUE
#define MATHICGB_UNREACHABLE
#define MATHICGB_RESTRICT

#endif

#include <cstddef>
#include <memory>

#ifdef MATHICGB_DEBUG
// don't force inline while debugging
#undef MATHICGB_INLINE
#define MATHICGB_INLINE inline

// we have to define DEBUG as lots of code assumes that asserts are turned
// on/off depending on DEBUG. Those should change to checking
// MATHICGB_DEBUG and then we can remove this define.
#define DEBUG
#include <iostream> // Useful for debugging.
#include <cassert>
#define MATHICGB_ASSERT(X) do{assert(X);}while(0)
#define MATHICGB_ASSERT_NO_ASSUME(X) MATHICGB_ASSERT(X)
#define MATHICGB_IF_DEBUG(X) X
#else
#define MATHICGB_ASSERT(X) MATHICGB_ASSUME(X)
#define MATHICGB_ASSERT_NO_ASSUME(X)
#define MATHICGB_IF_DEBUG(X)
#endif

#ifdef MATHICGB_SLOW_DEBUG
// for asserts that take a long time.
#define MATHICGB_SLOW_ASSERT(X) MATHICGB_ASSERT(X)
#else
#define MATHICGB_SLOW_ASSERT(X)
#endif

#include <utility>
/*
See http://herbsutter.com/gotw/_102/ for a reason to have a
make_unique function. It's pretty easy to do, too:

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

Unfortunately, MSVC does not have variadic templates, so this turns
into the monstrosity of overloads below. At least they got the perfect
forwarding working, otherwise this would have required N^2 overloads
for N parameters! Add more overloads below if you need more
parameters.
*/

template<class T>
std::unique_ptr<T> make_unique() {
  return std::unique_ptr<T>(new T());
}
template<class T, class A1>
std::unique_ptr<T> make_unique(A1&& a1) {
  return std::unique_ptr<T>(new T(std::forward<A1>(a1)));
}
template<class T, class A1, class A2>
std::unique_ptr<T> make_unique(A1&& a1, A2&& a2) {
  return std::unique_ptr<T>(new T(std::forward<A1>(a1), std::forward<A2>(a2)));
}
template<class T, class A1, class A2, class A3>
std::unique_ptr<T> make_unique(A1&& a1, A2&& a2, A3&& a3) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2), std::forward<A3>(a3)));
}
template<class T, class A1, class A2, class A3, class A4>
  std::unique_ptr<T> make_unique(A1&& a1, A2&& a2, A3&& a3, A4&& a4) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2),
           std::forward<A3>(a3), std::forward<A4>(a4)));
}
template<class T, class A1, class A2, class A3, class A4, class A5>
  std::unique_ptr<T> make_unique(A1&& a1, A2&& a2, A3&& a3, A4&& a4, A5&& a5) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2),
           std::forward<A3>(a3), std::forward<A4>(a4),
           std::forward<A5>(a5)));
}
template<class T, class A1, class A2, class A3, class A4, class A5, class A6>
  std::unique_ptr<T> make_unique
    (A1&& a1, A2&& a2, A3&& a3, A4&& a4, A5&& a5, A6&& a6) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2),
           std::forward<A3>(a3), std::forward<A4>(a4),
           std::forward<A5>(a5), std::forward<A6>(a6)));
}
template<class T, class A1, class A2, class A3, class A4, class A5, class A6,
class A7>
  std::unique_ptr<T> make_unique
    (A1&& a1, A2&& a2, A3&& a3, A4&& a4, A5&& a5, A6&& a6, A7&& a7) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2),
           std::forward<A3>(a3), std::forward<A4>(a4),
           std::forward<A5>(a5), std::forward<A6>(a6),
           std::forward<A7>(a7)));
}
template<class T, class A1, class A2, class A3, class A4, class A5, class A6,
class A7, class A8>
  std::unique_ptr<T> make_unique
    (A1&& a1, A2&& a2, A3&& a3, A4&& a4, A5&& a5, A6&& a6, A7&& a7,
    A8&& a8) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2),
           std::forward<A3>(a3), std::forward<A4>(a4),
           std::forward<A5>(a5), std::forward<A6>(a6),
           std::forward<A7>(a7), std::forward<A8>(a8)));
}
template<class T, class A1, class A2, class A3, class A4, class A5, class A6,
class A7, class A8, class A9>
  std::unique_ptr<T> make_unique
    (A1&& a1, A2&& a2, A3&& a3, A4&& a4, A5&& a5, A6&& a6, A7&& a7,
    A8&& a8, A9&& a9) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2),
           std::forward<A3>(a3), std::forward<A4>(a4),
           std::forward<A5>(a5), std::forward<A6>(a6),
           std::forward<A7>(a7), std::forward<A8>(a8),
           std::forward<A9>(a9)));
}
template<class T, class A1, class A2, class A3, class A4, class A5, class A6,
class A7, class A8, class A9, class A10>
  std::unique_ptr<T> make_unique
    (A1&& a1, A2&& a2, A3&& a3, A4&& a4, A5&& a5, A6&& a6, A7&& a7,
    A8&& a8, A9&& a9, A10&& a10) {
  return std::unique_ptr<T>
    (new T(std::forward<A1>(a1), std::forward<A2>(a2),
           std::forward<A3>(a3), std::forward<A4>(a4),
           std::forward<A5>(a5), std::forward<A6>(a6),
           std::forward<A7>(a7), std::forward<A8>(a8),
           std::forward<A9>(a9), std::forward<A9>(a10)));
}

template<class T>
std::unique_ptr<T[]> make_unique_array(size_t count) {
  return std::unique_ptr<T[]>(new T[count]);
}

// TODO: These types should be defined in some way that actually
// checks that these bit counts are right like in a configure script.
typedef unsigned long long uint64;
typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char uint8;

typedef signed long long int64;
typedef signed int int32;
typedef signed short int16;
typedef signed char int8;

/// Bizarrely, OpenMP 2 only supports signed loop
/// variables. This defect is fixed in OpenMP 3. MSVC 2012 only supports
/// OpenMP 2. So signed loop indices are forced for loops that are
/// parallelized using OpenMP.
typedef signed long OMPIndex;

static const size_t BitsPerByte = 8;
static const size_t MemoryAlignment = sizeof(void*);

/// The higher the value the more detailed output about what the program
/// is doing.
extern int tracingLevel;
