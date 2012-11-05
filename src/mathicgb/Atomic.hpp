#ifndef MATHICGB_ATOMIC_GUARD
#define MATHICGB_ATOMIC_GUARD

// We need this include for std::memory_order even if we are not
// using std::atomic.
#include <atomic>

namespace AtomicInternal {
  /// Class for deciding which implementation of atomic to use. The default is
  /// to use std::atomic which is a fine choice if std::atomic is implemented
  /// in a reasonable way by the standard library implementation you are using.
  template<class T, size_t size>
  struct ChooseAtomic {
    typedef std::atomic<T> type;
  };
}

#if defined(MATHICGB_USE_CUSTOM_ATOMIC_X86_X64_MSVC_4BYTE) || \
  defined(MATHICGB_USE_CUSTOM_ATOMIC_X86_X64_MSVC_8BYTE)
#include <Windows.h>
#undef max
#undef min
namespace AtomicInternal {
  /// Custom Atomic class. Use for sizes that are automatically atomic on
  /// the current platform.
  template<class T>
  class CustomAtomic {
  public:
    MATHICGB_INLINE CustomAtomic(): mValue() {}
    MATHICGB_INLINE CustomAtomic(T value): mValue(value) {}

    MATHICGB_INLINE T load(std::memory_order order) const {
      switch (order) {
      case std::memory_order_relaxed:
        // In this case there are no constraints on ordering, so all we need
        // to ensure is to perform an atomic read. That is already automatic
        // on aligned variables of this type.
        return mValue;

      case std::memory_order_consume: {
        // This order specifies to not move dependent reads above this load,
        // that is to we must not load *p before loading p.
        // The only mainstream CPU that will break this guarantee is DEC Alpha,
        // in which case this code should not be used. The compiler can also
        // speculate the value of p and load *p and then later check that p is
        // what it thought it was. That will break the guarantee that we need
        // to give, so we need a compiler optimization barrier that will
        // disable such speculative dependent read optimizations.
        //
        // Unfortunately on e.g. MSVC there is not a specific barrier for this
        // purpose so we end up blocking all read-movement, including
        // non-dependent reads. 
        const auto value = mValue;
        _ReadBarrier(); // MSVC
        return value;
      }

      case std::memory_order_acquire: {
        // This order specifies to not move any read above this load. Some CPUs
        // and all compilers can break this guarantee due to optimizations.
        // x86, ARM, SPARC and x64 has this guarantee automatically, though see
       //http://preshing.com/20120913/acquire-and-release-semantics#comment-20810
        // for an exception on x64 that I am going to ignore - if you use these
        // techniques it is up to you to ensure proper fencing. On other
        // platforms we need a hardware fence in addition to the compiler
        // optimization fence.
        const auto value = mValue;
        _ReadBarrier(); // MSVC
        return mValue;
      }

      case std::memory_order_seq_cst:
        // There must be some global order in which all sequentially consistent
        // atomic operations are considered to have happened in. This is automatic
        // on x64, ARM, SPARC and x64 too for reads (but not writes) - see:
        //   http://www.stdthread.co.uk/forum/index.php?topic=72.0
        _ReadBarrier(); // MSVC
        return mValue;

      case std::memory_order_release: // not available for load
      case std::memory_order_acq_rel: // not available for load
      default:
        MATHICGB_UNREACHABLE;
      }
    }

    MATHICGB_INLINE void store(const T value, std::memory_order order) {
      switch (order) {
      case std::memory_order_relaxed:
        // No ordering constraints and we are already assuming that loads
        // and stores to mValue are automatic so we can just store directly.
        mValue = value;
        break;

      case std::memory_order_release:
        // This ordering specifies that no prior writes are moved after
        // this write.
        _WriteBarrier();
        mValue = value;
        break;

      case std::memory_order_acq_rel:
        // This ordering specifies combined acquire and release guarantees:
        //  - no prior writes are moved after this operation
        //  - no later reads are moved before this operation
        // This requires CPU fencing on x86/x64 because normally reads can
        // be moved before writes to a different memory location by the CPU.
        _WriteBarrier();
        mValue = value;
        MemoryBarrier();
        break;

      case std::memory_order_seq_cst:
        // All operations happen in a globally consistent linear order.
        _WriteBarrier();
        mValue = value;
        MemoryBarrier();
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

#ifdef MATHICGB_USE_CUSTOM_ATOMIC_X86_X64_MSVC_4BYTE
  template<class T>
  struct ChooseAtomic<T, 4> {
    typedef CustomAtomic<T> type;
  };
#endif

#ifdef MATHICGB_USE_CUSTOM_ATOMIC_X86_X64_MSVC_8BYTE
  template<class T>
  struct ChooseAtomic<T, 8> {
    typedef CustomAtomic<T> type;
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
  MATHICGB_INLINE Atomic(): mValue() {}
  MATHICGB_INLINE Atomic(T value): mValue(value) {}

  MATHICGB_INLINE T load(
    std::memory_order order = std::memory_order_seq_cst
  ) const {
    MATHICGB_ASSERT(debugAligned());
    return mValue.load(order);
  }

  MATHICGB_INLINE void store(
    T value,
    std::memory_order order = std::memory_order_seq_cst
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
