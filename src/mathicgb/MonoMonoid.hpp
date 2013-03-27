#ifndef MATHICGB_MONO_MONOID_GUARD
#define MATHICGB_MONO_MONOID_GUARD

#include <vector>
#include <algorithm>
#include <memtailor.h>

/// Implements the monoid of (monic) monomials with integer
/// non-negative exponents. T must be an unsigned integer type that is
/// used to store each exponent of a monomial.
///
/// TODO: support grading and comparison.
template<class E>
class MonoMonoid {
public:
  // *** Types

  // Integer index representing a variable. Indices start at 0 and go
  // up to varCount() - 1 where varCount() is the number of variables.
  typedef size_t VarIndex;

  /// The type of each exponent of a monomial.
  typedef E Exponent;

  /// Iterator for the exponents in a monomial.
  typedef const Exponent* const_iterator;

  /// Represents a monomial and manages the memory underlying it. To
  /// refer to a non-owned monomial or to refer to a Mono, use MonoRef
  /// or ConstMonoRef. Do not use Mono& or Mono* if you do not have
  /// to, since that implies a double indirection when accessing the
  /// monomial.
  class Mono;

  /// A reference to a non-const monomial. Cannot be null, cannot be
  /// reassigned to refer to a different monomial and does not connote
  /// ownership - the same semantics as C++ references.
  class MonoRef;

  /// A reference to a monomial. As MonoRef, but you cannot change the
  /// monomial through this reference. Prefer this class over the
  /// other reference/pointer classes unless there is a reason not to.
  class ConstMonoRef;

  /// A pointer to a non-const monomial. Can be null and can be
  /// reassigned to refer to a different monomial - the same semantics
  /// as C++ pointers. Does not connote ownership.
  class MonoPtr;

  /// A pointer to a monomial. As MonoPtr, but you cannot change the
  /// monomial through this pointer.
  class ConstMonoPtr;

  /// A pool of memory for monomials.
  ///
  /// @todo: This approach is a poor fit for variable-sized
  /// monomials. So prefer other solutions where reasonable.
  class MonoPool;

  /// A vector of monomials. The interface is a subset of
  /// std::vector. Monomials can be appended (push_back). Only the
  /// last monomial can be mutated and monomials cannot be reordered
  /// or removed. These restrictions should make it easier to support
  /// variable-sized monomials in future. Change it if you need to
  /// break these restrictions, but first try to find an alternative.
  class MonoVector;


  // *** Constructors and accessors

  MonoMonoid(const VarIndex varCount): mVarCount(varCount) {}

  bool operator==(const MonoMonoid& monoid) const {
    return this == &monoid;
  }

  bool operator!=(const MonoMonoid& monoid) const {
    return !(*this == monoid);
  }

  /// Returns the number of variables. This is also the number of
  /// exponents in the exponent vector of a monomial.
  VarIndex varCount() const {return mVarCount;}


  // *** Monomial accessors and queries

  /// Returns iterator to the first exponent.
  const_iterator begin(ConstMonoRef mono) const {return mono.rawPtr();}

  /// Returns iterator to one-past-the-end of the range of exponents.
  const_iterator end(ConstMonoRef mono) const {
    return mono.rawPtr() + entriesPerMono();
  }

  /// Returns the exponent of var in mono.
  Exponent exponent(ConstMonoRef mono, const VarIndex var) const {
    MATHICGB_ASSERT(var < varCount());
    return mono.rawPtr()[var];
  } 

  bool equal(ConstMonoRef a, ConstMonoRef b) const {
    return std::equal(begin(a), end(a), begin(b));
  }

  /// Returns true if all the exponents of mono are zero. In other
  /// words, returns true if mono is the identity for multiplication
  /// of monomials.
  bool isIdentity(ConstMonoRef mono) const {
    return std::all_of(begin(mono), end(mono), [](Exponent e) {return e == 0;});
  }


  // *** Monomial mutating computations

  void copy(ConstMonoRef from, MonoRef to) const {
    std::copy_n(from.rawPtr(), entriesPerMono(), to.rawPtr());
  }

  void setExponent(const VarIndex var, const Exponent e, MonoRef mono) const {
    MATHICGB_ASSERT(var < varCount());
    mono.rawPtr()[var] = e;
  }

  void setIdentity(MonoRef mono) const {
    std::fill_n(mono.rawPtr(), entriesPerMono(), static_cast<Exponent>(0));
    MATHICGB_ASSERT(isIdentity(mono));
  }


  // *** Classes for holding and referring to monomials

  class ConstMonoPtr {
  public:
    ConstMonoPtr(): mMono(0) {}
    ConstMonoPtr(const ConstMonoPtr& mono): mMono(mono.rawPtr()) {}

    ConstMonoPtr operator=(const ConstMonoPtr& mono) {
      mMono = mono.mMono;
      return *this;
    }

    ConstMonoRef operator*() const {return *this;}

    bool isNull() const {return mMono == 0;}
    void toNull() {mMono = 0;}

  private:
    friend class MonoVector;
    friend class MonoPtr;
    friend class ConstMonoRef;

    const Exponent* rawPtr() const {return mMono;}
    ConstMonoPtr(const Exponent* mono): mMono(mono) {}

    const Exponent* mMono;
  };

  class MonoPtr {
  public:
    MonoPtr(): mMono(0) {}
    MonoPtr(const MonoPtr& mono): mMono(mono.rawPtr()) {}

    MonoPtr operator=(const MonoPtr& mono) {
      mMono = mono.mMono;
      return *this;
    }

    MonoRef operator*() const {return *this;}

    bool isNull() const {return mMono == 0;}
    void toNull() {mMono = 0;}

    operator ConstMonoPtr() const {return ConstMonoPtr(mMono);}

  private:
    friend class Mono;
    friend class MonoRef;
    friend class MonoVector;
    friend class MonoPool;

    Exponent* rawPtr() const {return mMono;}
    MonoPtr(Exponent* mono): mMono(mono) {}

    Exponent* mMono;
  };

  class Mono {
  public:
    Mono(): mMono(), mPool(0) {}

    Mono(Mono&& mono): mMono(mono.mMono), mPool(mono.mPool) {
      mono.mMono.toNull();
      mono.mPool = 0;
    }

    ~Mono() {toNull();}

    void operator=(Mono&& mono) {
      toNull();

      mMono = mono.mMono;
      mono.mMono.toNull();

      mPool = mono.mPool;
      mono.mPool = 0;
    }
    
    bool isNull() const {return mMono.isNull();}
    void toNull() {mPool->free(*this);}

    MonoPtr ptr() const {return mMono;}

    operator MonoRef() const {
      MATHICGB_ASSERT(!isNull());
      return *mMono;
    }

  private:
    Mono(const Mono&); // not available
    void operator=(const Mono&); // not available
    friend class MonoPool;

    Mono(const MonoPtr mono, MonoPool& pool):
      mMono(mono), mPool(&pool) {}

    Exponent* rawPtr() const {return mMono.rawPtr();}

    MonoPtr mMono;
    MonoPool* mPool;
  };

  class MonoRef {
  public:
    MonoPtr ptr() const {return mMono;}

    operator ConstMonoRef() const {return *static_cast<ConstMonoPtr>(mMono);}

  private:
    void operator=(const MonoRef&); // not available
    friend class MonoMonoid;
    friend class MonoPtr;

    MonoRef(MonoPtr mono): mMono(mono) {}
    Exponent* rawPtr() const {return mMono.rawPtr();}

    const MonoPtr mMono;
  };

  class ConstMonoRef {
  public:
    ConstMonoRef(const Mono& mono): mMono(mono.ptr()) {
      MATHICGB_ASSERT(!mono.isNull());
    }

    ConstMonoPtr ptr() const {return mMono;}

  private:
    void operator=(const MonoRef&); // not available
    friend class MonoMonoid;
    friend class ConstMonoPtr;

    ConstMonoRef(ConstMonoPtr mono): mMono(mono) {}
    const Exponent* rawPtr() const {return mMono.rawPtr();}

    const ConstMonoPtr mMono;
  };


  // *** Classes that provide memory resources for monomials

  class MonoPool {
  public:
    MonoPool(const MonoMonoid& monoid):
      mMonoid(monoid),
      mPool(sizeof(Exponent) * mMonoid.entriesPerMono()) {}

    Mono alloc() {
      const auto ptr = static_cast<Exponent*>(mPool.alloc());
      Mono mono(ptr, *this);
      monoid().setIdentity(mono);
      return mono;
    }

    void free(Mono& mono) {free(std::move(mono));}
    void free(Mono&& mono) {
      if (mono.isNull())
	return;
      mPool.free(mono.rawPtr());
      mono.mMono = 0;
      mono.mPool = 0;
    }

    const MonoMonoid& monoid() const {return mMonoid;}

  private:
    const MonoMonoid& mMonoid;
    memt::BufferPool mPool;
  };

  class MonoVector {
  private:
    typedef std::vector<Exponent> RawVector;

  public:
    /// Class for iterating through the monomials in a MonoVector.
    ///
    /// There is no operator->() since MonoRef does not have any
    /// relevant methods to call. Implement it if you need it.
    ///
    /// There are no postfix increment operator as prefix is
    /// better. Add it if you y need it (you probably do not).
    ///
    /// We could make this a random access iterator, but that would
    /// make it tricky to support variable-sized exponent vectors
    /// (e.g. sparse) in future and so far we have not needed random
    /// access.
    class const_iterator {
    public:
      typedef std::forward_iterator_tag iterator_category;
      typedef ConstMonoPtr value_type;
    
      const_iterator(): mIt(), mEntriesPerMono(0) {}
      const_iterator(const const_iterator& it):
        mIt(it.mIt), mEntriesPerMono(it.mEntriesPerMono) {}
    
      bool operator==(const const_iterator& it) const {return mIt == it.mIt;}
      bool operator!=(const const_iterator& it) const {return mIt != it.mIt;}

      ConstMonoRef operator*() {
	MATHICGB_ASSERT(debugValid());
	return *ConstMonoPtr(&*mIt);
      }

      const_iterator operator++() {
	MATHICGB_ASSERT(debugValid());
	mIt += mEntriesPerMono;
	return *this;
      }

    private:
      friend class MonoVector;
      bool debugValid() {return mEntriesPerMono > 0;}

      const_iterator(
        typename RawVector::const_iterator it,
	size_t entriesPerMono
      ): mIt(it), mEntriesPerMono(entriesPerMono) {}
      
      typename RawVector::const_iterator mIt;
      size_t mEntriesPerMono;		     
    };

    // ** Constructors and assignment
    MonoVector(const MonoMonoid& monoid): mMonoid(monoid) {}
    MonoVector(const MonoVector& v): mMonos(v.mMonos), mMonoid(v.monoid()) {}
    MonoVector(MonoVector&& v):
      mMonos(std::move(v.mMonos)), mMonoid(v.monoid()) {}

    MonoVector& operator=(const MonoVector& v) {
      MATHICGB_ASSERT(monoid() == v.monoid());
      mMonos = v.mMonos;
      return *this;
    }

    MonoVector& operator=(MonoVector&& v) {
      MATHICGB_ASSERT(monoid() == v.monoid());
      mMonos = std::move(v.mMonos);
      return *this;      
    }

    // ** Iterators
    const_iterator begin() const {
      return const_iterator(mMonos.begin(), mMonoid.entriesPerMono());
    }

    const_iterator end() const {
      return const_iterator(mMonos.end(), mMonoid.entriesPerMono());
    }

    const_iterator cbegin() const {return begin();}
    const_iterator cend() const {return end();}

    // ** Capacity
    size_t size() const {return mMonos.size() / monoid().entriesPerMono();}
    bool empty() const {return mMonos.empty();}

    // ** Element access
    ConstMonoRef front() const {
      MATHICGB_ASSERT(!empty());
      return *begin();
    }

    MonoRef back() {
      MATHICGB_ASSERT(!empty());
      const auto offset = mMonos.size() - monoid().entriesPerMono();
      return *MonoPtr(mMonos.data() + offset);
    }

    ConstMonoRef back() const {
      MATHICGB_ASSERT(!empty());
      const auto offset = mMonos.size() - monoid().entriesPerMono();
      return *ConstMonoPtr(mMonos.data() + offset);
    }

    // ** Modifiers
    void push_back(ConstMonoRef mono) {
      const auto offset = mMonos.size();
      mMonos.resize(offset + monoid().entriesPerMono());
      monoid().copy(mono, *MonoPtr(mMonos.data() + offset));
    }

    /// Appends the identity.
    void push_back() {
      const auto offset = mMonos.size();
      mMonos.resize(offset + monoid().entriesPerMono());      
    }

    void swap(MonoVector& v) {
      MATHICGB_ASSERT(&monoid() == &v.monoid());
      mMonos.swap(v.mMonos);
    }

    void clear() {mMonos.clear();}

    // ** Relational operators
    bool operator==(const MonoVector& v) const {
      MATHICGB_ASSERT(monoid() == v.monoid());
      return mMonos == v.mMonos;
    }
    bool operator!=(const MonoVector& v) const {return !(*this == v);}

    // ** Other
    size_t memoryBytesUsed() const {
      return mMonos.capacity() * sizeof(mMonos[0]);
    }

    const MonoMonoid& monoid() const {return mMonoid;}

  private:
    RawVector mMonos;
    const MonoMonoid& mMonoid;
  };


private:
  size_t entriesPerMono() const {return varCount();}

  const VarIndex mVarCount;
};

#endif
