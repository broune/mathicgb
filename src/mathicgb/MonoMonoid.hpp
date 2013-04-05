#ifndef MATHICGB_MONO_MONOID_GUARD
#define MATHICGB_MONO_MONOID_GUARD

#include <vector>
#include <algorithm>
#include <memtailor.h>
#include <type_traits>
#include <istream>
#include <ostream>
#include <cstdlib>
#include <cstring>
#include <mathic.h>

/// Implements the monoid of (monic) monomials with integer
/// non-negative exponents. Exponent must be an unsigned integer type that is
/// used to store each exponent of a monomial.
template<
  class Exponent,
  bool HasComponent = true,
  bool StoreHash = true,
  bool StoreOrder = true
>
class MonoMonoid;


template<class E, bool HC, bool SH, bool SO>
class MonoMonoid {
public:
  static_assert(std::numeric_limits<E>::is_signed, "");


  // *** Types

  // Integer index representing a variable. Indices start at 0 and go
  // up to varCount() - 1 where varCount() is the number of variables.
  typedef size_t VarIndex;

  /// The type of each exponent of a monomial.
  typedef E Exponent;

  /// Is true if the monomials come from a module.
  static const bool HasComponent = HC;

  /// Is true if the hash value is stored rather than computed at each 
  /// hash request. This imposes extra computation when updating a monomial,
  /// but for most operations that overhead is much less than the time for
  /// computing a hash value from scratch.
  static const bool StoreHash = SH;

  /// Is true if data to compare monomials is stored rather than computed
  /// at each comparison. As storeHash, there is overhead for this, but it
  /// is not much for most operations.
  static const bool StoreOrder = SO;

  /// Type used to indicate the component of a module monomial. For example,
  /// the component of xe_3 is 3.
  typedef typename std::make_unsigned<E>::type Component;

  /// Type used to store hash values of monomials.
  typedef typename std::make_unsigned<E>::type HashValue;

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

  // *** Temporary compatibility code for migrating off PolyRing
  friend class PolyRing;
  static MonoRef toRef(Exponent* e) {return MonoRef(e);}
  static ConstMonoRef toRef(const Exponent* e) {return ConstMonoRef(e);}
  static Exponent* toOld(MonoRef e) {return rawPtr(e);}
  static const Exponent* toOld(ConstMonoRef e) {return rawPtr(e);}

  

  // *** Constructors and accessors

  MonoMonoid(const std::vector<Exponent>& grading):
    mVarCount(grading.size()),
    mOrderEntryCount(StoreOrder),
    mOrderIndexBegin(HasComponent + mVarCount),
    mOrderIndexEnd(HasComponent + mVarCount + StoreOrder),
    mHashCoefficients(mVarCount),
    mGradingIsTotalDegree(
      [&]() {
        for (auto it = grading.begin(); it != grading.end(); ++it)
          if (*it != 1)
            return false;
        return true;
      }()
    ),
    mGrading(mGradingIsTotalDegree ? std::vector<Exponent>() : grading),
    mPool(*this)
  {
    if (!mGradingIsTotalDegree) {
      // Take negative values since reverse lex makes bigger values
      // give smaller monomials, but we need bigger degrees to give
      // bigger monomials.
      mGrading.resize(varCount());
      for (VarIndex var = 0; var < varCount(); ++var)
        mGrading[var] = -grading[var];
    }

    std::srand(0); // To use the same hash coefficients every time.
    for (VarIndex var = 0; var < varCount(); ++var)
      mHashCoefficients[var] = static_cast<HashValue>(std::rand());      
  }


  MonoMonoid(const VarIndex varCount):
    mVarCount(varCount),
    mOrderEntryCount(StoreOrder),
    mOrderIndexBegin(HasComponent + mVarCount),
    mOrderIndexEnd(HasComponent + mVarCount + StoreOrder),
    mHashCoefficients(mVarCount),
    mGradingIsTotalDegree(true),
    mGrading(),
    mPool(*this)
  {
    std::srand(0); // To use the same hash coefficients every time.
    for (VarIndex var = 0; var < varCount; ++var)
      mHashCoefficients[var] = static_cast<HashValue>(std::rand());      
  }

  /// Copies the ordering and varCount from monoid.
  template<class E2, bool HC2, bool SH2, bool SO2>
  MonoMonoid(const MonoMonoid<E2, HC2, SH2, SO2>& monoid):
    mVarCount(monoid.varCount()),
    mOrderEntryCount(StoreOrder),
    mOrderIndexBegin(HasComponent + mVarCount),
    mOrderIndexEnd(HasComponent + mVarCount + StoreOrder),
    mHashCoefficients(mVarCount),
    mGradingIsTotalDegree(monoid.mGradingIsTotalDegree),
    mGrading(monoid.mGrading),
    mPool(*this)
  {}
  // the compiler-generated copy constructor is a better match for overload
  // resolution than the template constructor above, so we need to provide
  // our own copy constructor to prevent the compiler-generated one from
  // being used.
  MonoMonoid(const MonoMonoid& monoid):
    mVarCount(monoid.varCount()),
    mOrderEntryCount(StoreOrder),
    mOrderIndexBegin(HasComponent + mVarCount),
    mOrderIndexEnd(HasComponent + mVarCount + StoreOrder),
    mHashCoefficients(mVarCount),
    mGradingIsTotalDegree(monoid.mGradingIsTotalDegree),
    mGrading(monoid.mGrading),
    mPool(*this)
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
    return ptr(mono, exponentsIndexBegin());
  }

  /// Returns iterator to one-past-the-end of the range of exponents.
  const_iterator end(ConstMonoRef mono) const {
    return ptr(mono, exponentsIndexEnd());
  }

  /// Returns the exponent of var in mono.
  Exponent exponent(ConstMonoRef mono, const VarIndex var) const {
    MATHICGB_ASSERT(var < varCount());
    return access(mono, exponentsIndexBegin() + var);
  } 

  /// Returns the component of the monomial. Monomials not from a
  /// module have component zero. In a module mono*e_i has component
  /// i. @todo: Have different monoids for module monomials and
  /// monomials and only offer this method for the module monomials.
  Component component(ConstMonoRef mono) const {
    MATHICGB_ASSERT(HasComponent);
    return access(mono, componentIndex());
  }

  /// Returns a hash value for the monomial. These are not guaranteed
  /// to be unique.
  HashValue hash(ConstMonoRef mono) const {
    MATHICGB_ASSERT(debugHashValid(mono));
    if (StoreHash)
      return static_cast<HashValue>(access(mono, hashIndex()));
    else
      return computeHash(mono);
  }

  /// Returns true if a and b are equal. Includes check for component.
  bool equal(ConstMonoRef a, ConstMonoRef b) const {
    for (auto i = entriesIndexBegin(); i != exponentsIndexEnd(); ++i)
      if (access(a, i) != access(b, i))
        return false;
    return true;
  }

  template<class MonoidA>
  bool equal(
    const MonoidA& monoidA,
    typename MonoidA::ConstMonoRef a,
    ConstMonoRef b
  ) const {
    // todo: assert compatible
    for (VarIndex var = 0; var < varCount(); ++var)
      if (monoidA.exponent(a, var) != exponent(b, var))
        return false;
    return true;
  }

  /// As equal(), but optimized for the case where true is returned.
  bool equalHintTrue(ConstMonoRef a, ConstMonoRef b) const {
    // if a[i] != b[i] then a[i] ^ b[i] != 0, so the or of all xors is zero
    // if and only if a equals b. This way we avoid having a branch to check
    // equality for every iteration of the loop, which is a win in the case
    // that none of the early-exit branches are taken - that is, when a equals
    // b.
    Exponent orOfXor = 0;
    for (VarIndex i = lastExponentIndex(); i != beforeEntriesIndexBegin(); --i)
      orOfXor |= access(a, i) ^ access(b, i);
    MATHICGB_ASSERT((orOfXor == 0) == equal(a, b));
    return orOfXor == 0;
  }

  bool isProductOf(
    ConstMonoRef a,
    ConstMonoRef b,
    ConstMonoRef ab
  ) const {
    for (VarIndex i = entriesIndexBegin(); i != exponentsIndexEnd(); ++i)
      if (access(ab, i) != access(a, i) + access(b, i))
        return false;
    return true;
  }

  bool isProductOfHintTrue(
    ConstMonoRef a, 
    ConstMonoRef b, 
    ConstMonoRef ab
  ) const {
    // We compare more than one exponent at a time using 64 bit integers. This 
    // might go one 32 bit value at the end too far, but since that space is
    // either a degree or a hash value that is fine --- those values will also
    // match if the monomials are equal. This does not work for negative
    // exponents since the overflowing bit will go into the next word.
    // It is OK that the degree field can be negative (a field we might go
    // into without caring about it because it shares a 64 bit field with
    // the last exponent), because it is at the end so the overflowing
    // bit will not interfere. For this reason we need to have a degree
    // or a hash value stored there - otherwise two equal monomials could
    // have different things stored next to them which would confuse this code.
    
    // todo: ensure 8 byte alignment. Though there seem to be no ill effects
    // for unaligned access. Performance seems to be no worse than for using
    // 32 bit integers directly.

    if (sizeof(Exponent) != 4 || (!StoreHash && !StoreOrder))
      return isProductOf(a, b, ab);

    uint64 orOfXor = 0;
    for (VarIndex i = varCount() / 2; i != beforeEntriesIndexBegin(); --i) {
      MATHICGB_ASSERT(access(a, i*2) >= 0);
      MATHICGB_ASSERT(i == varCount() / 2 || access(a, i*2+1) >= 0);
      
      uint64 A, B, AB;
      // We have to use std::memcpy here because just casting to a int64 breaks
      // the strict aliasing rule which implies undefined behavior. Both MSVC and
      // gcc don't actually call memcpy here. MSVC is a tiny bit slower for this
      // code than for casting while GCC seems to be exactly the same speed.
      std::memcpy(&A, ptr(a, i*2), 8);
      std::memcpy(&B, ptr(b, i*2), 8);
      std::memcpy(&AB, ptr(ab, i*2), 8);
      orOfXor |= AB ^ (A + B);
    }
    MATHICGB_ASSERT((orOfXor == 0) == isProductOf(a, b, ab));
    return orOfXor == 0; 
  }

  MATHICGB_INLINE bool isTwoProductsOfHintTrue(
    ConstMonoRef a1,
    ConstMonoRef a2,
    ConstMonoRef b,
    ConstMonoRef a1b,
    ConstMonoRef a2b
  ) const {
    if (sizeof(Exponent) != 4 || (!StoreHash && !StoreOrder))
      return (isProductOf(a1, b, a1b) && isProductOf(a2, b, a2b));

    uint64 orOfXor = 0;
    for (VarIndex i = varCount() / 2; i != beforeEntriesIndexBegin(); --i) {
      uint64 A1, A2, B, A1B, A2B;
      std::memcpy(&A1, ptr(a1, i*2), 8);
      std::memcpy(&A2, ptr(a2, i*2), 8);
      std::memcpy(&B, ptr(b, i*2), 8);
      std::memcpy(&A1B, ptr(a1b, i*2), 8);
      std::memcpy(&A2B, ptr(a2b, i*2), 8);
      orOfXor |= (A1B ^ (A1 + B)) | (A2B ^ (A2 + B));
    }
    MATHICGB_ASSERT
      ((orOfXor == 0) == (isProductOf(a1, b, a1b) && isProductOf(a2, b, a2b)));
    return orOfXor == 0;
  }

  /// Returns the hash of the product of a and b.
  HashValue hashOfProduct(ConstMonoRef a, ConstMonoRef b) const {
    // See computeHash() for an explanation of all the casts.
    const auto hashA = static_cast<HashValue>(hash(a));
    const auto hashB = static_cast<HashValue>(hash(b));
    return static_cast<HashValue>(static_cast<Exponent>(hashA + hashB));
  }

  /// Returns true if all the exponents of mono are zero. In other
  /// words, returns true if mono is the identity for multiplication
  /// of monomials.
  bool isIdentity(ConstMonoRef mono) const {
    return std::all_of(begin(mono), end(mono), [](Exponent e) {return e == 0;});
  }

  /// Returns true if a divides b. Equal monomials divide each other.
  bool divides(ConstMonoRef div, ConstMonoRef into) const {
    // todo: enable this when the code works with it
    //if (HasComponent && component(div) != component(into))
    //  return false;
    for (auto i = exponentsIndexBegin(); i < exponentsIndexEnd(); ++i)
      if (access(div, i) > access(into, i))
        return false;
    return true;
  }

  template<class MonoidA>
  bool divides(
    const MonoidA& monoidA,
    typename MonoidA::ConstMonoRef a,
    ConstMonoRef b
  ) const {
    // todo: fix other divisibility functions to work properly for component too.
    MATHICGB_ASSERT(monoidA.varCount() == varCount());
    MATHICGB_ASSERT(!MonoidA::HasComponent || HasComponent);
    MATHICGB_ASSERT(monoidA.debugValid(a));
    MATHICGB_ASSERT(debugValid(b));
    // todo: enable this when the code works with it
    //if (HasComponent && component(div) != component(into))
    //  return false;
    //if (
    //  MonoidA::HasComponent &&
    //  HasComponent &&
    //  monoidA.component(a) != component(b)
    //)
    //  return false;
    for (VarIndex var = 0; var < varCount(); ++var)
      if (monoidA.exponent(a, var) > exponent(b, var))
        return false;
    return true;
  }

  /// Returns true if div divides lcm(a, b).
  bool dividesLcm(ConstMonoRef div, ConstMonoRef a, ConstMonoRef b) const {
    MATHICGB_ASSERT(debugLcmCheck(*this, a, *this, b));
    MATHICGB_ASSERT(debugValid(div));

    for (auto i = exponentsIndexBegin(); i != exponentsIndexEnd(); ++i) {
      const auto dive = access(div, i);
      if (access(div, i) > access(a, i) && access(div, i) > access(b, i))
        return false;
    }
    return true;
  }

  template<class MonoidDiv, class MonoidA>
  bool dividesLcm(
    const MonoidDiv& monoidDiv,
    typename MonoidDiv::ConstMonoRef div,
    const MonoidA& monoidA,
    typename MonoidA::ConstMonoRef a,
    ConstMonoRef b
  ) const {
    MATHICGB_ASSERT(monoidDiv.debugLcmCheck(monoidA, a, *this, b));
    MATHICGB_ASSERT(monoidDiv.debugValid(div));

    for (VarIndex var = 0; var < varCount(); ++var) {
      const auto e = monoidDiv.exponent(div, var);
      if (e > monoidA.exponent(a, var) && e > exponent(b, var))
        return false;
    }
    return true;
  }

  /// Returns true if lcm(a,b) == lcmAB.
  bool isLcm(ConstMonoRef a, ConstMonoRef b, ConstMonoRef lcmAB) const {
    MATHICGB_ASSERT(debugLcmCheck(*this, a, *this, b));
    MATHICGB_ASSERT(debugValid(lcmAB));

    for (auto i = exponentsIndexBegin(); i != exponentsIndexEnd(); ++i)
      if (access(lcmAB, i) != std::max(access(a, i), access(b, i)))
        return false;
    return true;
  }

  template<class MonoidA, class MonoidB>
  bool isLcm(
    const MonoidA& monoidA,
    typename MonoidA::ConstMonoRef a,
    const MonoidB& monoidB,
    typename MonoidB::ConstMonoRef b,
    ConstMonoRef lcmAB
  ) const {
    MATHICGB_ASSERT(debugLcmCheck(monoidA, a, monoidB, b));
    MATHICGB_ASSERT(debugValid(lcmAB));

    if (HasComponent) {
      if (MonoidA::HasComponent) {
        if (monoidA.component(a) != component(lcmAB))
          return false;
      } else {
        MATHICGB_ASSERT(MonoidB::HasComponent);
        if (monoidB.component(b) != component(lcmAB))
          return false;
      }
    }

    for (VarIndex var = 0; var < varCount(); ++var) {
      if (
        ptr(lcmAB, exponentsIndexBegin())[var] !=
        std::max(monoidA.exponent(a, var), monoidB.exponent(b, var))
      )
        return false;
    }
    return true;
  }

  // Graded reverse lexicographic order. The grading is total degree.
  CompareResult compare(ConstMonoRef a, ConstMonoRef b) const {
    MATHICGB_ASSERT(debugOrderValid(a));
    MATHICGB_ASSERT(debugOrderValid(b));

    const auto cmp = degree(a) - degree(b);
    if (cmp < 0) return GreaterThan;
    if (cmp > 0) return LessThan;

    for (auto i = exponentsIndexEnd() - 1; i != beforeEntriesIndexBegin(); --i) {
      const auto cmp = access(a, i) - access(b, i);
      if (cmp < 0) return GreaterThan;
      if (cmp > 0) return LessThan;
    }
    return EqualTo;
  }

  bool lessThan(ConstMonoRef a, ConstMonoRef b) const {
    return compare(a, b) == LessThan;
  }

  /// Returns true if gcd(a, b) == 1.
  bool relativelyPrime(ConstMonoRef a, ConstMonoRef b) const {
    for (auto i = exponentsIndexBegin(); i != exponentsIndexEnd(); ++i)
      if (access(a, i) > 0 && access(b, i) > 0)
        return false;
    return true;
  }

  // If this method returns true for monomials a and b then it is
  // guaranteed that multiplying a and b together will not overflow
  // the integers in the representation.
  bool hasAmpleCapacity(ConstMonoRef mono) const {
    const auto halfMin = std::numeric_limits<Exponent>::min() / 2;
    const auto halfMax = std::numeric_limits<Exponent>::max() / 2;
    MATHICGB_ASSERT(halfMin <= 0);
    const auto limit = std::min(-halfMin, halfMax);
    const auto inRange = [&](Exponent value)
      {return -limit <= value && value <= limit;};

    for (VarIndex i = exponentsIndexBegin(); i != exponentsIndexEnd(); ++i)
      if (!inRange(access(mono, i)))
        return false;
    if (!inRange(degree(mono)))
      return false;
    return true;
  }

  /// Returns the degree of mono using the grading on the monoid.
  Exponent degree(ConstMonoRef mono) const {
    MATHICGB_ASSERT(debugOrderValid(mono));
    if (StoreOrder)
      return access(mono, orderIndexBegin());
    else
      return computeDegree(mono);
  }


  // *** Monomial mutating computations

  /// Copes the parameter from to the parameter to.
  void copy(ConstMonoRef from, MonoRef to) const {
    MATHICGB_ASSERT(debugValid(from));

    std::copy_n(rawPtr(from), entryCount(), rawPtr(to));

    MATHICGB_ASSERT(debugValid(to));
  }

  template<class MonoidFrom>
  void copy(
    const MonoidFrom& monoidFrom,
    typename MonoidFrom::ConstMonoRef from,
    MonoRef to
  ) const {
    // todo: extract this in checker method
    MATHICGB_ASSERT(HasComponent == MonoidFrom::HasComponent);
    MATHICGB_ASSERT(monoidFrom.debugValid(from));
    MATHICGB_ASSERT(monoidFrom.varCount() == varCount());
    MATHICGB_ASSERT
      ((std::is_same<Exponent, typename MonoidFrom::Exponent>::value));

    if (HasComponent)
      access(to, componentIndex()) = monoidFrom.component(from);
    for (VarIndex var = 0; var < varCount(); ++var)
      ptr(to, exponentsIndexBegin())[var] = monoidFrom.exponent(from, var);
    if (StoreOrder)
      access(to, orderIndexBegin()) = monoidFrom.degree(from);
    if (StoreHash)
      access(to, hashIndex()) = monoidFrom.hash(from);

    MATHICGB_ASSERT(debugValid(to));
    // todo: check equal
  }

  /// Set the exponent of var to newExponent in mono.
  void setExponent(
    const VarIndex var,
    const Exponent newExponent,
    MonoRef mono
  ) const {
    MATHICGB_ASSERT(var < varCount());

    auto& exponent = access(mono, exponentsIndexBegin() + var);
    const auto oldExponent = exponent;
    exponent = newExponent;

    updateOrderData(var, oldExponent, newExponent, mono);
    updateHashExponent(var, oldExponent, newExponent, mono);

    MATHICGB_ASSERT(debugValid(mono));
  }

  /// Sets all the exponents of mono. exponents must point to an array
  /// of size varCount().
  void setExponents(const Exponent* exponents, MonoRef mono) const {
    MATHICGB_ASSERT(exponents != 0);

    if (HasComponent)
      access(mono, componentIndex()) = 0;
    std::copy_n(exponents, varCount(), ptr(mono, exponentsIndexBegin()));
    setOrderData(mono);
    setHash(mono);

    MATHICGB_ASSERT(debugValid(mono));
  }

  /// Sets mono to 1, which is the identity for multiplication.
  void setIdentity(MonoRef mono) const {
    std::fill_n(rawPtr(mono), entryCount(), static_cast<Exponent>(0));

    MATHICGB_ASSERT(debugValid(mono));
    MATHICGB_ASSERT(isIdentity(mono));
  }

  /// Sets the component of mono to newComponent.
  void setComponent(Component newComponent, MonoRef mono) const {
    MATHICGB_ASSERT(HasComponent);

    auto& component = access(mono, componentIndex());
    const auto oldComponent = component;
    component = newComponent;
    updateHashComponent(oldComponent, newComponent, mono);

    MATHICGB_ASSERT(debugValid(mono));
  }

  /// Sets prod to a*b.
  void multiply(ConstMonoRef a, ConstMonoRef b, MonoRef prod) const {
    MATHICGB_ASSERT(debugValid(a));
    MATHICGB_ASSERT(debugValid(b));

    for (auto i = lastEntryIndex(); i != beforeEntriesIndexBegin(); --i)
      access(prod, i) = access(a, i) + access(b, i);

    MATHICGB_ASSERT(debugValid(prod));
  }

  /// Sets prod to a*prod.
  void multiplyInPlace(ConstMonoRef a, MonoRef prod) const {
    MATHICGB_ASSERT(debugValid(a));
    MATHICGB_ASSERT(debugValid(prod));

    for (auto i = entriesIndexBegin(); i < entriesIndexEnd(); ++i)
      access(prod, i) += access(a, i);

    MATHICGB_ASSERT(debugValid(prod));      
  }

  /// Sets quo to num/by. by must divide num.
  void divide(ConstMonoRef by, ConstMonoRef num, MonoRef quo) const {
    MATHICGB_ASSERT(divides(by, num));
    MATHICGB_ASSERT(debugValid(num));
    MATHICGB_ASSERT(debugValid(by));

    for (auto i = entriesIndexBegin(); i < entriesIndexEnd(); ++i)
      access(quo, i) = access(num, i) - access(by, i);

    MATHICGB_ASSERT(debugValid(quo));
  }

  /// Sets num to num/by. by must divide num.
  void divideInPlace(ConstMonoRef by, MonoRef num) const {
    MATHICGB_ASSERT(divides(by, num));
    MATHICGB_ASSERT(debugValid(by));
    MATHICGB_ASSERT(debugValid(num));

    for (auto i = entriesIndexBegin(); i < entriesIndexEnd(); ++i)
      access(num, i) -= access(by, i);

    MATHICGB_ASSERT(debugValid(num));
  }

  /// Sets quo to num/by. If by does not divide num then quo will have
  /// negative exponents.
  void divideToNegative(ConstMonoRef by, ConstMonoRef num, MonoRef quo) const {
    MATHICGB_ASSERT(debugValid(num));
    MATHICGB_ASSERT(debugValid(by));
    MATHICGB_ASSERT(
      !HasComponent ||
      component(by) == 0 ||
      component(by) == component(num)
    );

    for (auto i = entriesIndexBegin(); i < entriesIndexEnd(); ++i)
      access(quo, i) = access(num, i) - access(by, i);

    MATHICGB_ASSERT(debugValid(quo));
  }

  /// Sets aColonB to a:b and bColonA to b:a.
  void colons(
    ConstMonoRef a,
    ConstMonoRef b,
    MonoRef aColonB,
    MonoRef bColonA
  ) const {
    MATHICGB_ASSERT(debugValid(a));
    MATHICGB_ASSERT(debugValid(b));

    if (HasComponent) {
      MATHICGB_ASSERT(component(a) == component(b));
      access(aColonB, componentIndex()) = 0;
      access(bColonA, componentIndex()) = 0;
    }

    for (auto i = exponentsIndexBegin(); i != exponentsIndexEnd(); ++i) {
      const auto ae = access(a, i);
      const auto be = access(b, i);
      const auto max = std::max(ae, be);
      access(aColonB, i) = max - be;
      access(bColonA, i) = max - ae;
    }
    setOrderData(aColonB);
    setHash(aColonB);
    setOrderData(bColonA);
    setHash(bColonA);

    MATHICGB_ASSERT(debugValid(aColonB));
    MATHICGB_ASSERT(debugValid(bColonA));
  }

  /// Sets lcmAB to the lcm of a and b.
  void lcm(ConstMonoRef a, ConstMonoRef b, MonoRef lcmAB) const {
    if (HasComponent) {
      MATHICGB_ASSERT(component(a) == component(b));
      access(lcmAB, componentIndex()) = access(a, componentIndex());
    }
    for (auto i = exponentsIndexBegin(); i != exponentsIndexEnd(); ++i)
      access(lcmAB, i) = std::max(access(a, i), access(b, i));
    setOrderData(lcmAB);
    setHash(lcmAB);

    MATHICGB_ASSERT(debugValid(lcmAB));
    MATHICGB_ASSERT(isLcm(a, b, lcmAB));
  }

  template<class MonoidA, class MonoidB>
  void lcm(
    const MonoidA& monoidA,
    typename MonoidA::ConstMonoRef a,
    const MonoidB& monoidB,
    typename MonoidB::ConstMonoRef b,
    MonoRef lcmAB
  ) const {
    MATHICGB_ASSERT(debugLcmCheck(monoidA, a, monoidB, b));

    if (HasComponent) {
      access(lcmAB, componentIndex()) =
        MonoidA::HasComponent ? monoidA.component(a) : monoidB.component(b);
    }

    for (VarIndex var = 0; var < varCount(); ++var) {
      ptr(lcmAB, exponentsIndexBegin())[var] =
        std::max(monoidA.exponent(a, var), monoidB.exponent(b, var));
    }

    setOrderData(lcmAB);
    setHash(lcmAB);

    MATHICGB_ASSERT(debugValid(lcmAB));
    MATHICGB_ASSERT(isLcm(monoidA, a, monoidB, b, lcmAB));
  }

  Mono alloc() const {return mPool.alloc();}
  void free(Mono&& mono) const {mPool.free(std::move(mono));}

  /// Parses a monomial out of a string. Valid examples: 1 abc a2bc
  /// aA. Variable names are case sensitive. Whitespace terminates the
  /// parse as does any other character that is not a letter or a
  /// digit.  The monomial must not include a coefficient, not even 1,
  /// except for the special case of a 1 on its own. An input like 1a
  /// will be parsed as two separate monomials. A suffix like <2> puts
  /// the monomial in component 2, so a5<2> is a^5e_2. The default
  /// component is 0.
  void parseM2(std::istream& in, MonoRef mono) const;

  // Inverse of parseM2().
  void printM2(ConstMonoRef mono, std::ostream& out) const;


  // *** Classes for holding and referring to monomials

  class ConstMonoPtr {
  public:
    ConstMonoPtr(): mMono(0) {}
    ConstMonoPtr(const ConstMonoPtr& mono): mMono(rawPtr(mono)) {}

    ConstMonoPtr operator=(const ConstMonoPtr& mono) {
      mMono = mono.mMono;
      return *this;
    }

    ConstMonoRef operator*() const {return *this;}

    bool isNull() const {return mMono == 0;}
    void toNull() {mMono = 0;}

  private:
    friend class MonoMonoid;

    const Exponent* internalRawPtr() const {return mMono;}
    ConstMonoPtr(const Exponent* mono): mMono(mono) {}

    const Exponent* mMono;
  };

  class MonoPtr {
  public:
    MonoPtr(): mMono(0) {}
    MonoPtr(const MonoPtr& mono): mMono(mono.mMono) {}

    MonoPtr operator=(const MonoPtr& mono) {
      mMono = mono.mMono;
      return *this;
    }

    MonoRef operator*() const {return *this;}

    bool isNull() const {return mMono == 0;}
    void toNull() {mMono = 0;}

    operator ConstMonoPtr() const {return ConstMonoPtr(mMono);}

  private:
    friend class MonoMonoid;

    Exponent* internalRawPtr() const {return mMono;}
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
    friend class MonoMonoid;

    Mono(const MonoPtr mono, MonoPool& pool):
      mMono(mono), mPool(&pool) {}

    Exponent* internalRawPtr() const {return rawPtr(mMono);}

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

    MonoRef(MonoPtr mono): mMono(mono) {}
    Exponent* internalRawPtr() const {return rawPtr(mMono);}

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

    ConstMonoRef(ConstMonoPtr mono): mMono(mono) {}
    const Exponent* internalRawPtr() const {return rawPtr(mMono);}

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
      mPool.free(rawPtr(mono));
      mono.mMono = 0;
      mono.mPool = 0;
    }

    const MonoMonoid& monoid() const {return mMonoid;}

  private:
    MonoPool(const MonoPool&); // not available
    void operator=(const MonoPool&); // not available

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
      MATHICGB_ASSERT(monoid().debugValid(mono));
      const auto offset = mMonos.size();
      mMonos.resize(offset + monoid().entryCount());
      monoid().copy(mono, *MonoPtr(mMonos.data() + offset));
    }

    /// Appends the identity.
    void push_back() {
      const auto offset = mMonos.size();
      mMonos.resize(offset + monoid().entryCount());
      MATHICGB_ASSERT(monoid().isIdentity(back()));
      MATHICGB_ASSERT(monoid().debugValid(back()));
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
  // Grants access to other template instantiations.
  template<class E2, bool HC2, bool SH2, bool SO2>
  friend class MonoMonoid;

  // The main point here is to grant access to rawPtr().
  friend class Mono;
  friend class MonoRef;
  friend class ConstMonoRef;
  friend class MonoPtr;
  friend class ConstMonoPtr;
  friend class MonoVector;
  friend class MonoPool;

  bool debugValid(ConstMonoRef mono) const {
    MATHICGB_ASSERT(debugOrderValid(mono));
    MATHICGB_ASSERT(debugHashValid(mono));
    return true;
  }

  template<class MonoidA, class MonoidB>
  bool debugLcmCheck(
    const MonoidA& monoidA,
    typename MonoidA::ConstMonoRef a,
    const MonoidB& monoidB,
    typename MonoidB::ConstMonoRef b
  ) const {
    MATHICGB_ASSERT(monoidA.varCount() == varCount());
    MATHICGB_ASSERT(monoidB.varCount() == varCount());
    MATHICGB_ASSERT
      ((std::is_same<Exponent, typename MonoidA::Exponent>::value));
    MATHICGB_ASSERT
      ((std::is_same<Exponent, typename MonoidB::Exponent>::value));
    MATHICGB_ASSERT
      (HasComponent == (MonoidA::HasComponent || MonoidB::HasComponent));
    MATHICGB_ASSERT(monoidA.debugValid(a));
    MATHICGB_ASSERT(monoidB.debugValid(b));
    MATHICGB_ASSERT(
      !HasComponent ||
      !MonoidA::HasComponent ||
      !MonoidB::HasComponent ||
      monoidA.component(a) == monoidB.component(b)
  );

    return true;
  }

  // *** Accessing fields of a monomial
  template<class M>
  static auto rawPtr(M&& m) -> decltype(m.internalRawPtr()) {
    return m.internalRawPtr();
  }

  Exponent* ptr(MonoRef& m, const VarIndex index) const {
    MATHICGB_ASSERT(index <= entryCount());
    return rawPtr(m) + index;
  }

  const Exponent* ptr(ConstMonoRef& m, const VarIndex index) const {
    MATHICGB_ASSERT(index <= entryCount());
    return rawPtr(m) + index;
  }

  Exponent& access(MonoRef& m, const VarIndex index) const {
    MATHICGB_ASSERT(index < entryCount());
    return rawPtr(m)[index];
  }

  const Exponent& access(ConstMonoRef& m, const VarIndex index) const {
    MATHICGB_ASSERT(index < entryCount());
    return rawPtr(m)[index];
  }

  // *** Implementation of monomial ordering
  bool debugOrderValid(ConstMonoRef mono) const {
#ifdef MATHICGB_DEBUG
    if (!StoreOrder)
      return true;
    // Check assumptions for layout in memory.
    MATHICGB_ASSERT(orderIndexBegin() == exponentsIndexEnd());
    MATHICGB_ASSERT(orderIndexBegin() + orderEntryCount() == orderIndexEnd());
    MATHICGB_ASSERT(orderIndexEnd() <= entryCount());
    MATHICGB_ASSERT(orderEntryCount() == 1);

    if (mGradingIsTotalDegree)
      MATHICGB_ASSERT(mGrading.empty());
    else
      MATHICGB_ASSERT(mGrading.size() == varCount());

    // Check the order data of mono
    MATHICGB_ASSERT
      (rawPtr(mono)[orderIndexBegin()] == computeDegree(mono));
#endif
    return true;
  }

  void setOrderData(MonoRef mono) const {
    if (!StoreOrder)
      return;

    rawPtr(mono)[orderIndexBegin()] = computeDegree(mono);      
    MATHICGB_ASSERT(debugOrderValid(mono));
  }

  void updateOrderData(
    const VarIndex var,
    const Exponent oldExponent,
    const Exponent newExponent,
    MonoRef mono
  ) const {
    if (!StoreOrder)
      return;

    MATHICGB_ASSERT(var < varCount());
    if (mGradingIsTotalDegree)
      rawPtr(mono)[orderIndexBegin()] += oldExponent - newExponent;
    else {
      MATHICGB_ASSERT(mGrading.size() == varCount());
      rawPtr(mono)[orderIndexBegin()] -=
        mGrading[var] * (oldExponent - newExponent);
    }
    MATHICGB_ASSERT(debugOrderValid(mono));
  }

  Exponent computeDegree(ConstMonoRef mono) const {
    Exponent degree = 0;
    const auto end = this->end(mono);
    if (mGradingIsTotalDegree) {
      for (auto it = begin(mono); it != end; ++it)
        degree -= *it;
    } else {
      for (auto var = 0; var < varCount(); ++var)
        degree += access(mono, exponentsIndexBegin() + var) * mGrading[var];
    }
    return degree;
  }


  // *** Implementation of hash value computation

  bool debugHashValid(ConstMonoRef mono) const {
    if (!StoreHash)
      return true;

    // We cannot call hash() here since it calls this method.
    MATHICGB_ASSERT(hashIndex() < entryCount());
    // todo: we cannot make this check right now because the legacy
    // integration with PolyRing can create monomials with unset hash.
    // MATHICGB_ASSERT(rawPtr(mono)[hashIndex()] == computeHash(mono));
    return true;
  }

  HashValue computeHash(ConstMonoRef mono) const {
    HashValue hash = HasComponent ? component(mono) : 0;
    for (VarIndex var = 0; var < varCount(); ++var)
      hash += static_cast<HashValue>(exponent(mono, var)) * mHashCoefficients[var];

    // Hash values are stored as exponents. If the cast to an exponent
    // changes the value, then we need computeHashValue to match that
    // change by casting to an exponent and back. Otherwise the computed
    // hash value will not match a hash value that has been stored.
    return static_cast<HashValue>(static_cast<Exponent>(hash));
  }

  void setHash(MonoRef mono) const {
    if (!StoreHash)
      return;
    rawPtr(mono)[hashIndex()] = computeHash(mono);
    MATHICGB_ASSERT(debugHashValid(mono));
  }

  void updateHashComponent(
    const Exponent oldComponent,
    const Exponent newComponent,
    MonoRef mono
  ) const {
    if (!StoreHash)
      return;
    rawPtr(mono)[hashIndex()] += newComponent - oldComponent;
    MATHICGB_ASSERT(debugHashValid(mono));
  }

  void updateHashExponent(
    const VarIndex var,
    const Exponent oldExponent,
    const Exponent newExponent,
    MonoRef mono
  ) const {
    if (!StoreHash)
      return;
    MATHICGB_ASSERT(var < varCount());
    rawPtr(mono)[hashIndex()] +=
      (newExponent - oldExponent) * mHashCoefficients[var];
    MATHICGB_ASSERT(debugHashValid(mono));
  }


  // *** Code determining the layout of monomials in memory
  // Layout in memory:
  //   [component] [exponents...] [order data...] [hash]

  /// Returns how many Exponents are necessary to store a
  /// monomial. This can include other data than the exponents, so
  /// this number can be larger than varCount().
  size_t entryCount() const {return mOrderIndexEnd + StoreHash;}

  /// Returns how many Exponents are necessary to store the extra data
  /// used to compare monomials quickly.
  size_t orderEntryCount() const {
    // static_assert(StoreOrder, "");
    return mOrderEntryCount;
  }

  VarIndex componentIndex() const {
    //static_assert(HasComponent, "");
    return 0;
  }

  VarIndex exponentsIndexBegin() const {return HasComponent;}
  VarIndex exponentsIndexEnd() const {return exponentsIndexBegin() + varCount();}
  VarIndex lastExponentIndex() const {return exponentsIndexEnd() - 1;}

  VarIndex orderIndexBegin() const {
    //static_assert(StoreOrder, "");
    return mOrderIndexBegin;
  }

  VarIndex orderIndexEnd() const {
    //static_assert(StoreOrder, "");
    return mOrderIndexEnd;
  }

  VarIndex hashIndex() const {
    //static_assert(StoreHash, "");
    return mOrderIndexEnd;
  }

  VarIndex entriesIndexBegin() const {return 0;}
  VarIndex entriesIndexEnd() const {return entryCount();}
  VarIndex beforeEntriesIndexBegin() const {return entriesIndexBegin() - 1;}
  VarIndex lastEntryIndex() const {return entriesIndexEnd() - 1;}

  const VarIndex mVarCount;
  const VarIndex mOrderEntryCount;
  const VarIndex mOrderIndexBegin;
  const VarIndex mOrderIndexEnd;

  /// Take dot product of exponents with this vector to get hash value.
  std::vector<HashValue> mHashCoefficients;

  /// This is initialized before mGrading, so it has to be ordered
  /// above mGrading.
  const bool mGradingIsTotalDegree;

  /// The degree of a monomial is the dot product of the exponent
  /// vector with this vector. If mGradingIsTotalDegree is true then
  /// mGrading is empty but implicitly it is a vector of ones.
  std::vector<Exponent> mGrading;

  mutable MonoPool mPool;
};

namespace MonoMonoidHelper {
  /// ostream and istream handle characters differently from other
  /// integers. Use unchar to cast chars to a different type that get
  /// handled as other integers do.
  template<class T>
  struct unchar {typedef int type;};

  // Yes: char, signed char and unsigned char are 3 distinct types.
  template<>
  struct unchar<char> {typedef short type;};
  template<>
  struct unchar<signed char> {typedef short type;};
  template<>
  struct unchar<unsigned char> {typedef unsigned short type;};
}

template<class E, bool HC, bool SH, bool SO>
void MonoMonoid<E, HC, SH, SO>::parseM2(std::istream& in, MonoRef mono) const {
  using MonoMonoidHelper::unchar;
  // todo: signal error on exponent overflow

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
    auto& exponent = access(mono, exponentsIndexBegin() + var);
    if (isdigit(in.peek())) {
      typename unchar<Exponent>::type e;
      in >> e;
      exponent = static_cast<Exponent>(e);
    } else
      exponent = 1;
    sawSome = true;
  }

  if (in.peek() == '<') {
    if (!HasComponent) {
      mathic::reportError("Read unexpected < for start of module component\n");
      return;
    }

    in.get();
    if (!isdigit(in.peek())) {
      mathic::reportError("Component was not integer.");
      return;
    }
    typename unchar<Exponent>::type e;
    in >> e;
    access(mono, componentIndex()) = static_cast<Exponent>(e);
    if (in.peek() != '>') {
      mathic::reportError("Component < was not matched by >.");
      return;
    }
    in.get();
  }

  setOrderData(mono);
  setHash(mono);
  MATHICGB_ASSERT(debugValid(mono));
}

template<class E, bool HC, bool SH, bool SO>
void MonoMonoid<E, HC, SH, SO>::printM2(
  ConstMonoRef mono,
  std::ostream& out
) const {
  using MonoMonoidHelper::unchar;
  const auto letterCount = 'z' - 'a' + 1;

  bool printedSome = false;
  for (VarIndex var = 0; var < varCount(); ++var) {
    if (exponent(mono, var) == 0)
      continue;
    char letter;
    if (var < letterCount)
      letter = 'a' + static_cast<char>(var);
    else if (var < 2 * letterCount)
      letter = 'A' + (static_cast<char>(var) - letterCount);
    else {
      mathic::reportError("Too few letters in alphabet to print variable.");
      return;
    }
    printedSome = true;
    out << letter;
    if (exponent(mono, var) != 1)
      out << static_cast<typename unchar<Exponent>::type>(exponent(mono, var));
  }
  if (!printedSome)
    out << '1';
  if (HasComponent && component(mono) != 0) {
    out << '<'
        << static_cast<typename unchar<Exponent>::type>(component(mono))
        << '>';
  }
}

#endif
