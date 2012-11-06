#ifndef MATHICGB_ATOMIC_GUARD
#define MATHICGB_ATOMIC_GUARD

// We need this include for std::memory_order even if we are not
// using std::atomic.
#include <atomic>

#if defined(_MSC_VER) && defined(MATHICGB_USE_CUSTOM_ATOMIC_X86_X64)

/// Tells the compiler (not the CPU) to not reorder reads across this line.
#define MATHICGB_COMPILER_READ_MEMORY_BARRIER _ReadBarrier()

/// Tells the compiler (not the CPU) to not reorder writes across this line.
#define MATHICGB_COMPILER_WRITE_MEMORY_BARRIER _WriteBarrier()

/// Tells the compiler (not the CPU) to not reorder reads and writes across
/// this line.
#define MATHICGB_COMPILER_READ_WRITE_MEMORY_BARRIER _ReadWriteBarrier()

/// Tells the CPU and also the compiler to not reorder reads and writes
/// across this line.
#define MATHICGB_CPU_READ_WRITE_MEMORY_BARRIER MemoryBarrier()

/// Loads a variable with sequentially consistent ordering. The variable
/// must be aligned and have a size such that aligned loads of that size
/// are atomic.
#define MATHICGB_SEQ_CST_LOAD(REF) \
  ::AtomicInternalMsvc::SeqCstSelect<decltype(REF)>::load(REF)

/// Stores a variable with sequentially consistent ordering. The variable
/// must be aligned and have a size such that aligned reads of that size
/// are atomic.
#define MATHICGB_SEQ_CST_STORE(VALUE, REF) \
  ::AtomicInternalMsvc::SeqCstSelect<decltype(REF)>::store(VALUE, REF)

#include <Windows.h>
// Windows.h defines macroes max and min that mess up things like std::max and
// std::numeric_limits<T>::max. So we need to undefine those macroes.
#undef max
#undef min
namespace AtomicInternalMsvc {
  template<class T, size_t size> struct SeqCst {};
#ifdef MATHICGB_USE_CUSTOM_ATOMIC_4BYTE
  template<class T> struct SeqCst<T, 4> {
    static T load(const T& ref) {
      return (T)_InterlockedOr((volatile LONG*)&ref, 0);
    }
    static void store(const T value, T& ref) {
      _InterlockedExchange((volatile LONG*)&ref, (LONG)value);
    }
  };
#endif
#ifdef MATHICGB_USE_CUSTOM_ATOMIC_8BYTE
  template<class T> struct SeqCst<T, 8> {
    static T load(const T& ref) {
      return (T)_InterlockedOr64((volatile _LONGLONG*)&ref, 0);
    }
    static void store(const T value, T& ref) {
      _InterlockedExchange64((volatile _LONGLONG*)&ref, (_LONGLONG)value);
    }
  };
#endif
  template<class T> struct SeqCstSelect : public SeqCst<T, sizeof(T)> {};
}
#endif

#if defined(__GNUC__) && defined(MATHICGB_USE_CUSTOM_ATOMIC_X86_X64)

// As far as I can tell this is not documented to work, but it is the
// only way to do this on GCC and it is what the Linux kernel does, so
// that will have to be good enough for me.
#define MATHICGB_COMPILER_READ_WRITE_MEMORY_BARRIER \
  __asm__ __volatile__ ("" ::: "memory");

// As far as I can tell there is no way to do a partial optimization
// barrier on GCC, so we have to do the full barrier every time.
#define MATHICGB_COMPILER_READ_MEMORY_BARRIER \
  MATHICGB_COMPILER_READ_WRITE_MEMORY_BARRIER

#define MATHICGB_COMPILER_WRITE_MEMORY_BARRIER \
  MATHICGB_COMPILER_READ_WRITE_MEMORY_BARRIER

#define MATHICGB_CPU_READ_WRITE_MEMORY_BARRIER __sync_synchronize()

#define MATHICGB_SEQ_CST_LOAD(REF) \
  AtomicInternalGCC::SeqCst<decltype(REF)>::load(REF)
#define MATHICGB_SEQ_CST_STORE(VALUE, REF) \
  AtomicInternalGCC::SeqCst<decltype(REF)>::store(VALUE, REF)

namespace AtomicInternalGCC {
  template<class T> struct SeqCst {
    static T load(const T& ref) {
      const auto ptr = static_cast<volatile T*>(const_cast<T*>(&ref));
      return __sync_fetch_and_or((volatile T*)&ref, 0);
    }
    static void store(const T value, T& ref) {
      const auto ptr = static_cast<volatile T*>(&ref);
      while (!__sync_bool_compare_and_swap(ptr, *ptr, value)) {}
    }
  };
  template<class T> struct SeqCst<const T> : public SeqCst<T> {};
}

#endif

namespace AtomicInternal {
#ifdef MATHICGB_USE_FAKE_ATOMIC
  // This class has the same interface as the actual custom atomic
  // class but it does absolutely no synchronization and it does not
  // constrain compiler optimizations in any way. The purpose of this class
  // is to enable it while running single core to compare the single core
  // overhead of the atomic ordering constraints.
  template<class T>
  class FakeAtomic {
  public:
    FakeAtomic(): mValue() {}
    FakeAtomic(T value): mValue(value) {}
    T load(const std::memory_order) const {return mValue;}
    void store(const T value, const std::memory_order order) {mValue = value;}

  private:
    T mValue;
  };

  template<class T, size_t size>
  struct ChooseAtomic {
    typedef FakeAtomic<T> type;
  };

#else
  /// Class for deciding which implementation of atomic to use. The default is
  /// to use std::atomic which is a fine choice if std::atomic is implemented
  /// in a reasonable way by the standard library implementation you are using.
  template<class T, size_t size>
  struct ChooseAtomic {
    typedef std::atomic<T> type;
  };
#endif
}

#ifdef MATHICGB_USE_CUSTOM_ATOMIC_X86_X64
namespace AtomicInternal {
  /// Custom Atomic class for x86 and x64. Uses special compiler instructions
  /// for barriers. Only instantiate this for sizes where aligned reads and
  /// writes are guaranteed to be atomic - this class only takes care of the
  /// ordering constraints using CPU and compiler fences. Since the directives
  /// to achieve this are coming from the compiler it is very strange that
  /// any compiler ships with a std::atomic that is worse than this - but
  /// that is very much the case.
  ///
  /// There are 5 kinds of reorderings that we are concerned with here. Let
  /// S,S' be stores and let L,L' be stores. Note that these short-hands may
  /// be idiosyncratic - feel free to find some standard terminology from
  /// some prominent source and fix this to reflect that.
  ///
  ///   SS: Store-after-store: Reorder S,S' to S',S
  ///   SL: Store-after-load: Reorder S,L to L,S
  ///   LS: Load-after-store: Reorder L,S to S,L
  ///   LL: Load-after-load: Reorder L,L' to L',L
  ///   DLL: Dependent-load-after-load: As LL but L' depends on L. For example
  ///     reordering the load of p->a to before the load of p is a DLL.
  ///
  /// The DEC Alpha processor will perform all of these reorderings in the
  /// absense of memory barriers telling it not to do that, including DLL.
  /// DLL can happen on DEC Alpha if p->a is cached locally while p is not.
  /// Then p will be loaded from memory while p->a is loaded from the cache,
  /// which is functionally identical to loading p->a before p since we may
  /// see a value of p->a that was stored before the value of p. This happens
  /// even if the processor that stored p did a full memory barrier between
  /// storing p->a and storing p.
  ///
  /// Compilers will also perform all of these reorderings to optimize the
  /// code - even including DLL. DLL happens if the compiler guesses what
  /// the value of p is, loads p->a and then checks that the guess for p
  /// was correct. This directly causes p->a to be actually loaded before p.
  /// These kinds of optimizations turn up in profile-driven optimization,
  /// but it is always allowed unless we tell the compiler not to do it.
  ///
  /// You can check this out here:
  ///   http://en.wikipedia.org/wiki/Memory_ordering
  ///
  /// On x86 and x64 only SL is doe by the CPU, so we need a CPU barrier to
  /// prevent that and nothing else. The compiler is free to perform all of
  /// these reorderings, so we need lots of compiler optimization barriers
  /// to deal with all of these cases.
  ///
  /// Some of the quotes below are from
  ///
  ///   http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1525.htm
  template<class T>
  class CustomAtomicX86X64 {
  public:
    CustomAtomicX86X64(): mValue() {}
    CustomAtomicX86X64(T value): mValue(value) {}

    MATHICGB_INLINE
    T load(const std::memory_order order) const {
      switch (order) {
      case std::memory_order_relaxed:
        // The only constraint here is that if you read *p, then you will never
        // after that read a value of *p that was stored before the value
        // you just read, where "before" is in terms of either the same thread
        // that did the writing or external synchronization of another thread
        // with this thread. This is automaticaly guaranteed on this platform
        // and the compiler cannot break this guarantee.
        return mValue;

      case std::memory_order_consume: {
        // Loads in this thread that depend on the loaded value must not be
        // reordered to before this load. So no DLL reorderings past this
        // load from after to before (up). So we need a read barrier AFTER the
        // load. It is a compiler only barrier since the CPU does not do DLL
        // reorderings. 
        const auto value = mValue;
        MATHICGB_COMPILER_READ_MEMORY_BARRIER;
        return value;
      }

      case std::memory_order_acquire: {
        // Loads in this thread must not be reordered to before this load.
        // So no LL reorderings past this load from after to before (up).
        // So we need a barrier AFTER the load. It is a compiler only barrier
        // since the CPU does not do LL reorderings.
        const auto value = mValue;
        MATHICGB_COMPILER_READ_MEMORY_BARRIER;
        return mValue;
      }

      case std::memory_order_seq_cst:
        // There must be some global order in which all sequentially consistent
        // atomic operations are considered to have happened in. This is automatic
        // on x64, ARM, SPARC and x64 too for reads (but not writes) - see:
        //   http://www.stdthread.co.uk/forum/index.php?topic=72.0
        return MATHICGB_SEQ_CST_LOAD(mValue);

      case std::memory_order_release: // not available for load
      case std::memory_order_acq_rel: // not available for load
      default:
        MATHICGB_UNREACHABLE;
      }
    }

    MATHICGB_INLINE
    void store(const T value, const std::memory_order order) {
      switch (order) {
      case std::memory_order_relaxed:
        // No ordering constraints here other than atomicity and as noted
        // for relaxed load so we can just store directly.
        mValue = value;
        break;

      case std::memory_order_release:
        // Stores in this thread must not be reordered to after this store.
        // So no SS reorderings past this load from before to after (down).
        // So we need a barrier BEFORE the load. It is a compiler only barrier
        // since the CPU does not do SS reorderings.
        MATHICGB_COMPILER_WRITE_MEMORY_BARRIER;
        mValue = value;
        break;

      case std::memory_order_acq_rel:
        // Combine the guarantees for std::memory_order_acquire and
        // std::memory_order_release. So no loads moved up past here (SL) and
        // no stores moved down past here (LL). We need a compiler barrier
        // BEFORE the load to avoid LL and a CPU barrier (implies also a
        // compiler barrier AFTER the load to avoid SL, since the CPU can in
        // fact do SL reordering.
        MATHICGB_COMPILER_WRITE_MEMORY_BARRIER;
        mValue = value;
        MATHICGB_CPU_READ_WRITE_MEMORY_BARRIER;
        break;

      case std::memory_order_seq_cst:
        // All operations happen in a globally consistent linear order. I am
        // sure if this can be achieved with barriers but I know that it can be
        // achieved with locked instructions, so I am using that.
        MATHICGB_SEQ_CST_STORE(value, mValue);
        break;

      case std::memory_order_consume: // not available for store
      case std::memory_order_acquire: // not available for store
      default:
        MATHICGB_UNREACHABLE;
      }
    }

  private:
    T mValue;
  };

#ifdef MATHICGB_USE_CUSTOM_ATOMIC_4BYTE
  template<class T>
  struct ChooseAtomic<T, 4> {
    typedef CustomAtomicX86X64<T> type;
  };
#endif

#ifdef MATHICGB_USE_CUSTOM_ATOMIC_8BYTE
  template<class T>
  struct ChooseAtomic<T, 8> {
    typedef CustomAtomicX86X64<T> type;
  };
#endif
}
#endif

/// This class is equivalent to std::atomic<T*>. Some functions from the
/// interface of std::atomic are missing - add them as necessary. Do not add
/// operator= and operator T() --- it is better to make the code explicit
/// about when and how loading and storing of atomic variables occurs.
///
/// The purpose of the class is that it performs far better than
/// std::atomic for some implementations. For example the std::atomic in MSVC
/// 2012 performs a compare-and-swap operation on a load even with the
/// paramter std::memory_order_relaxed.
///
/// We force all the functions to be inline because they can contain switches
/// on the value of std::memory_order. This will usually be a compile-time
/// constant parameter so that after inlining the switch will disappear. Yet
/// the code size of the switch may make some compilers avoid the inline.
template<class T>
class Atomic {
public:
  Atomic(): mValue() {}
  Atomic(T value): mValue(value) {}

  MATHICGB_INLINE
  T load(const std::memory_order order = std::memory_order_seq_cst) const {
    MATHICGB_ASSERT(debugAligned());
    return mValue.load(order);
  }

  MATHICGB_INLINE
  void store(
    const T value,
    const std::memory_order order = std::memory_order_seq_cst
  ) {
    MATHICGB_ASSERT(debugAligned());
    mValue.store(value, order);
  }

private:
  Atomic(const Atomic<T>&); // not available
  void operator=(const Atomic<T>&); // not available

  bool debugAligned() const {
    return reinterpret_cast<size_t>(&mValue) % sizeof(void*) == 0;
  }

  typename AtomicInternal::ChooseAtomic<T, sizeof(T)>::type mValue;
};

#endif
