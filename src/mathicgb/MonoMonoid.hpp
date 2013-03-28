#ifndef MATHICGB_MONO_MONOID_GUARD
#define MATHICGB_MONO_MONOID_GUARD

#include <vector>
#include <algorithm>
#include <memtailor.h>
#include <type_traits>
#include <istream>
#include <ostream>
#include <mathic.h>

/// Implements the monoid of (monic) monomials with integer
/// non-negative exponents. T must be an unsigned integer type that is
/// used to store each exponent of a monomial.
///
/// TODO: support grading and comparison.
template<class E>
class MonoMonoid {
public:
  static_assert(std::numeric_limits<E>::is_signed, "");


  // *** Types

  // Integer index representing a variable. Indices start at 0 and go
  // up to varCount() - 1 where varCount() is the number of variables.
  typedef size_t VarIndex;

  /// The type of each exponent of a monomial.
  typedef E Exponent;

  /// Type used to indicate the component of a module monomial. For example,
  /// the component of xe_3 is 3.
  typedef typename std::make_unsigned<E>::type Component;

  /// Type used to store hash values of monomials.
  //typedef typename std::make_unsigned<E>::type HashValue;

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

  /// For indicating the result of comparing one monomial to another.
  enum CompareResult {
    LessThan = -1,
    EqualTo = 0,
    GreaterThan = 1
  };


  // *** Constructors and accessors

  MonoMonoid(const VarIndex varCount):
    mVarCount(varCount),
    mOrderEntryCount(1),
    mOrderIndexBegin(1 + varCount),
    mOrderIndexEnd(2 + varCount)
  {}

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
  const_iterator begin(ConstMonoRef mono) const {
    return mono.rawPtr() + exponentsIndexBegin();
  }

  /// Returns iterator to one-past-the-end of the range of exponents.
  const_iterator end(ConstMonoRef mono) const {
    return mono.rawPtr() + exponentsIndexEnd();
  }

  /// Returns the exponent of var in mono.
  Exponent exponent(ConstMonoRef mono, const VarIndex var) const {
    MATHICGB_ASSERT(var < varCount());
    return begin(mono)[var];
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

  Component component(ConstMonoRef mono) const {
    return *(mono.rawPtr() + componentIndex());
  }

  // Graded reverse lexicographic order. The grading is total degree.
  CompareResult compare(ConstMonoRef a, ConstMonoRef b) const {
    MATHICGB_ASSERT(debugOrderValid(a));
    MATHICGB_ASSERT(debugOrderValid(b));

    const auto rawA = rawPtr(a);
    const auto rawB = rawPtr(b);
    const auto stop = exponentsIndexBegin() - 1;
    for (auto var = orderIndexEnd() - 1; var != stop; --var) {
      const auto cmp = rawA[var] - rawB[var];
      if (cmp < 0) return GreaterThan;
      if (cmp > 0) return LessThan;
    }
    return EqualTo;
  }

  bool lessThan(ConstMonoRef a, ConstMonoRef b) const {
    return compare(a, b) == LessThan;
  }


  // *** Monomial mutating computations

  void copy(ConstMonoRef from, MonoRef to) const {
    std::copy_n(from.rawPtr(), entryCount(), to.rawPtr());
  }

  void setExponent(const VarIndex var, const Exponent e, MonoRef mono) const {
    MATHICGB_ASSERT(var < varCount());
    auto& exponent = (rawPtr(mono) + exponentsIndexBegin())[var];
    const auto oldExponent = exponent;
    exponent = e;
    updateOrderInformation(var, oldExponent, exponent, mono);
  }

  void setIdentity(MonoRef mono) const {
    std::fill_n(mono.rawPtr(), entryCount(), static_cast<Exponent>(0));
    MATHICGB_ASSERT(isIdentity(mono));
  }

  void setComponent(Component comp, MonoRef mono) const {
    *(mono.rawPtr() + componentIndex()) = comp;
  }

  /// Parses a monomial out of a string. Valid examples: 1 abc a2bc
  /// aA. Variable names are case sensitive. Whitespace terminates the
  /// parse as does any other character that is not a letter or a
  /// digit.  The monomial must not include a coefficient, not even 1,
  /// except for the special case of a 1 on its own. An input like 1a
  /// will be parsed as two separate monomials. A suffix like <2> puts
  /// the monomial in component 2, so a5<2> is a^5e_2. The default
  /// component is 0.
  void parseM2(std::istream& in, MonoRef mono) const {
    setIdentity(mono);

    bool sawSome = false;
    while (true) {
      const char next = in.peek();
      if (!sawSome && next == '1') {
	in.get();
	break;
      }

      VarIndex var;
      const auto letterCount = 'z' - 'a' + 1;
      if ('a' <= next && next <= 'z')
	var = next - 'a';
      else if ('A' <= next && next <= 'Z')
	var = (next - 'A') + letterCount;
      else if (sawSome)
	break;
      else {
	mathic::reportError("Could not parse monomial.");
	return;
      }
      MATHICGB_ASSERT(var < 2 * letterCount);
      if (var >= varCount()) {
	mathic::reportError("Unknown variable.");
	return;
      }
      in.get();
      auto& exponent = rawPtr(mono)[exponentsIndexBegin() + var];
      if (isdigit(in.peek()))
	in >> exponent;
      else
	exponent = 1;
      sawSome = true;
    }

    if (in.peek() == '<') {
      in.get();
      if (!isdigit(in.peek())) {
	mathic::reportError("Component was not integer.");
	return;
      }
      in >> rawPtr(mono)[componentIndex()];
      if (in.peek() != '>') {
	mathic::reportError("Component < was not matched by >.");
	return;
      }
      in.get();
    }

    setOrderData(mono);
  }

  // Inverse of parseM2().
  void printM2(ConstMonoRef mono, std::ostream& out) const {
    const auto letterCount = 'z' - 'a' + 1;

    bool printedSome = false;
    for (VarIndex var = 0; var < varCount(); ++var) {
      const auto exponent = begin(mono)[var];
      if (exponent == 0)
	continue;
      char letter;
      if (var < letterCount)
	letter = 'a' + var;
      else if (var < 2 * letterCount)
	letter = 'A' + (var - letterCount);
      else {
	mathic::reportError("Not enough letters in alphabet to print variable.");
	return;
      }
      printedSome = true;
      out << letter;
      if (exponent != 1)
	out << exponent;    
    }
    if (!printedSome)
      out << '1';
    if (component(mono) != 0)
      out << '<' << component(mono) << '>';
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
      mPool(sizeof(Exponent) * mMonoid.entryCount()) {}

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
	size_t entryCount
      ): mIt(it), mEntriesPerMono(entryCount) {}
      
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
      return const_iterator(mMonos.begin(), mMonoid.entryCount());
    }

    const_iterator end() const {
      return const_iterator(mMonos.end(), mMonoid.entryCount());
    }

    const_iterator cbegin() const {return begin();}
    const_iterator cend() const {return end();}

    // ** Capacity
    size_t size() const {return mMonos.size() / monoid().entryCount();}
    bool empty() const {return mMonos.empty();}

    // ** Element access
    ConstMonoRef front() const {
      MATHICGB_ASSERT(!empty());
      return *begin();
    }

    MonoRef back() {
      MATHICGB_ASSERT(!empty());
      const auto offset = mMonos.size() - monoid().entryCount();
      return *MonoPtr(mMonos.data() + offset);
    }

    ConstMonoRef back() const {
      MATHICGB_ASSERT(!empty());
      const auto offset = mMonos.size() - monoid().entryCount();
      return *ConstMonoPtr(mMonos.data() + offset);
    }

    // ** Modifiers
    void push_back(ConstMonoRef mono) {
      const auto offset = mMonos.size();
      mMonos.resize(offset + monoid().entryCount());
      monoid().copy(mono, *MonoPtr(mMonos.data() + offset));
    }

    /// Appends the identity.
    void push_back() {
      const auto offset = mMonos.size();
      mMonos.resize(offset + monoid().entryCount());      
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

    /// As parseM2 on monoid, but accepts a non-empty space-separated
    /// list of monomials. The monomials are appended to the end of
    /// the vector.
    void parseM2(std::istream& in) {
      while(true) {
	push_back();
	monoid().parseM2(in, back());
	if (in.peek() != ' ')
	  break;
	in.get();
      }
    }

    /// The inverse of parseM2.
    void printM2(std::ostream& out) const {
      for (auto it = begin(); it != end(); ++it) {
	if (it != begin())
	  out << ' ';
	monoid().printM2(*it, out);
      }
      out << '\n';
    }

    const MonoMonoid& monoid() const {return mMonoid;}

  private:
    RawVector mMonos;
    const MonoMonoid& mMonoid;
  };


private:
  template<class M>
  static auto rawPtr(M&& m) -> decltype(m.rawPtr()) {return m.rawPtr();}

  // *** Implementation of monomial ordering
  bool debugOrderValid(ConstMonoRef mono) const {
#ifdef MATHICGB_DEBUG
    // Check assumptions for layout in memory.
    MATHICGB_ASSERT(orderIndexBegin() == exponentsIndexEnd());
    MATHICGB_ASSERT(orderIndexBegin() + orderEntryCount() == orderIndexEnd());
    MATHICGB_ASSERT(orderEntryCount() == 1);

    // Check the order data of mono
    MATHICGB_ASSERT
      (rawPtr(mono)[orderIndexBegin()] == -computeTotalDegree(mono));
#endif
    return true;
  }

  void setOrderData(MonoRef mono) const {
    rawPtr(mono)[orderIndexBegin()] = -computeTotalDegree(mono);      
    MATHICGB_ASSERT(debugOrderValid(mono));
  }

  void updateOrderInformation(
    const VarIndex var,
    const Exponent oldExponent,
    const Exponent newExponent,
    MonoRef mono
  ) const {
    rawPtr(mono)[orderIndexBegin()] += oldExponent - newExponent;
    MATHICGB_ASSERT(debugOrderValid(mono));
  }


  Exponent computeTotalDegree(ConstMonoRef mono) const {
    Exponent degree = 0;
    const auto end = this->end(mono);
    for (auto it = begin(mono); it != end; ++it)
      degree += *it;
    return degree;
  }

  // *** Code describing the layout of monomials in memory
  // Layout in memory:
  //   [component] [exponents...] [order data...] [hash]

  /// Returns how many Exponents are necessary to store a
  /// monomial. This can include other information than the exponents,
  /// so this number can be larger than varCount().
  size_t entryCount() const {return mOrderIndexEnd;}

  /// Returns how many Exponents are necessary to store the extra data
  /// used to compare monomials quickly.
  size_t orderEntryCount() const {return mOrderEntryCount;}

  VarIndex componentIndex() const {return 0;}
  VarIndex exponentsIndexBegin() const {return 1;}
  VarIndex exponentsIndexEnd() const {return varCount() + 1;}
  VarIndex orderIndexBegin() const {return mOrderIndexBegin;}
  VarIndex orderIndexEnd() const {return mOrderIndexEnd;}
  //VarIndex hashIndex() const {return mOrderIndexEnd;

  const VarIndex mVarCount;
  const VarIndex mOrderEntryCount;
  const VarIndex mOrderIndexBegin;
  const VarIndex mOrderIndexEnd;

  /// Take dot product of exponents with this vector to get hash value.
  //std::vector<HashValue> mHashVals;

};

#endif
