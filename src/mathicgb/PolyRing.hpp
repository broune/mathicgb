// Copyright 2011 Michael E. Stillman

#ifndef _polyRing_h_
#define _polyRing_h_

#define MATHICGB_USE_MONOID
#ifdef MATHICGB_USE_MONOID
#include "MonoMonoid.hpp"

#define MATHICGB_USE_FIELD
#include "PrimeField.hpp"
#endif

#include <assert.h>
#include <string>
#include <vector>
#include <memtailor.h>
#include <cstdio>
#include <cstring>
#include <limits>


#define LT (-1)
#define EQ 0
#define GT 1

#ifdef MATHICGB_USE_FIELD
template<class T>
PrimeField<
  typename std::make_unsigned<
    typename std::remove_reference<T>::type
  >::type
> makeField(T charac) {
  return charac;
}
#endif

/** Returns a^-1 mod modulus. It is required that 0 < a < modulus. */
template<class T>
T modularInverse(T a, T modulus) {
#ifdef MATHICGB_USE_FIELD
  MATHICGB_ASSERT(0 < a);
  MATHICGB_ASSERT(a < modulus);
  auto f = makeField(modulus);
  return f.inverse(f.toElementInRange(a)).value();
#else
  // we do two turns of the extended Euclidian algorithm per
  // loop. Usually the sign of x changes each time through the loop,
  // but we avoid that by representing every other x as its negative,
  // which is the value minusLastX. This way no negative values show
  // up.
  MATHICGB_ASSERT(0 < a);
  MATHICGB_ASSERT(a < modulus);
#ifdef MATHICGB_DEBUG
  T origA = a;
#endif
  T b = modulus;
  T minusLastX = 0;
  T x = 1;
  while (true) {
    MATHICGB_ASSERT(x <= modulus);
    MATHICGB_ASSERT(minusLastX <= modulus);

    // first turn
    if (a == 1)
      break;
    const T firstQuotient = b / a;
    b -= firstQuotient * a;
    minusLastX += firstQuotient * x;

    // second turn
    if (b == 1) {
      MATHICGB_ASSERT(minusLastX != 0);
      MATHICGB_ASSERT(minusLastX < modulus);
      x = modulus - minusLastX;
      break;
    }
    const T secondQuotient = a / b;
    a -= secondQuotient * b;
    x += secondQuotient * minusLastX;
  }
  MATHICGB_ASSERT(x >= 1);
  MATHICGB_ASSERT(x < modulus);
  MATHICGB_ASSERT_NO_ASSUME((static_cast<uint64>(origA) * x) % modulus == 1);
  return x;
#endif
}

template<class T>
struct ModularProdType {};
template<> struct ModularProdType<uint8> {typedef uint16 type;};
template<> struct ModularProdType<uint16> {typedef uint32 type;};
template<> struct ModularProdType<uint32> {typedef uint64 type;};
template<> struct ModularProdType<int8> {typedef int16 type;};
template<> struct ModularProdType<int16> {typedef int32 type;};
template<> struct ModularProdType<int32> {typedef int64 type;};

/** Returns a*b mod modulus.  It is required that 0 <= a, b < modulus. */
template<class T>
T modularProduct(T a, T b, T modulus) {
#ifdef MATHICGB_USE_FIELD
  MATHICGB_ASSERT(0 < a);
  MATHICGB_ASSERT(a < modulus);
  MATHICGB_ASSERT(0 <= b);
  MATHICGB_ASSERT(b < modulus);
  auto f = makeField(modulus);
  return f.product(f.toElementInRange(a), f.toElementInRange(b)).value();
#else
  typedef typename ModularProdType<T>::type BigT;
  MATHICGB_ASSERT(0 <= a);
  MATHICGB_ASSERT(a < modulus);
  MATHICGB_ASSERT(0 <= b);
  MATHICGB_ASSERT(b < modulus);
  BigT bigProd = static_cast<BigT>(a) * b;
  MATHICGB_ASSERT(a == 0 || bigProd / a == b);
  return static_cast<T>(bigProd % modulus);
#endif
}

/** Returns a+b mod modulus.  It is required that 0 <= a, b < modulus. */
template<class T>
T modularSum(T a, T b, T modulus) {
#ifdef MATHICGB_USE_FIELD
  MATHICGB_ASSERT(0 < a);
  MATHICGB_ASSERT(a < modulus);
  MATHICGB_ASSERT(0 <= b);
  MATHICGB_ASSERT(b < modulus);
  auto f = makeField(modulus);
  return f.sum(f.toElementInRange(a), f.toElementInRange(b)).value();
#else
  typedef typename ModularProdType<T>::type BigT;
  MATHICGB_ASSERT(0 <= a);
  MATHICGB_ASSERT(a < modulus);
  MATHICGB_ASSERT(0 <= b);
  MATHICGB_ASSERT(b < modulus);
  BigT bigSum = static_cast<BigT>(a) + b;
  MATHICGB_ASSERT(bigSum - a == b);
  return static_cast<T>(bigSum % modulus);
#endif
} 

/** Returns -a mod modulus. It is required that 0 <= a < modulus. */
template<class T>
T modularNegative(T a, T modulus) {
#ifdef MATHICGB_USE_FIELD
  MATHICGB_ASSERT(0 < a);
  MATHICGB_ASSERT(a < modulus);
  auto f = makeField(modulus);
  return f.negative(f.toElementInRange(a)).value();
#else
  MATHICGB_ASSERT(0 <= a);
  MATHICGB_ASSERT(a < modulus);
  return a == 0 ? 0 : modulus - a;
#endif
}

/** Returns -a mod modulus. It is required that 0 < a < modulus. */
template<class T>
T modularNegativeNonZero(T a, T modulus) {
#ifdef MATHICGB_USE_FIELD
  MATHICGB_ASSERT(0 < a);
  MATHICGB_ASSERT(a < modulus);
  auto f = makeField(modulus);
  return f.negativeNonZero(f.toElementInRange(a)).value();
#else
  MATHICGB_ASSERT(0 < a);
  MATHICGB_ASSERT(a < modulus);
  return modulus - a;
#endif
}

typedef int32 exponent ;
typedef uint32 HashValue;
typedef long coefficient;
#ifdef MATHICGB_USE_MONOID
typedef MonoMonoid<exponent> Monoid;
#endif
#ifdef MATHICGB_USE_FIELD
typedef PrimeField<unsigned long> Field;
#endif

typedef exponent* vecmonomial; // includes a component
typedef coefficient const_coefficient;

class Monomial;

class ConstMonomial
{
  friend class PolyRing;
  friend class OrderA;
  friend class OrderB;
  friend class OrderC;
  friend class OrderD;
  friend class OrderE;
public:
  //* Wrap a raw pointer to create a monomial
  // Assumptions:
  //  1. This is done in the presence of a PolyRing
  //  2. Space for the monomial has been created
  ConstMonomial() : mValue(0) {}
  ConstMonomial(const exponent *val) : mValue(val) {}

  inline const Monomial& castAwayConst() const;

  bool isNull() const { return mValue == 0; }

  exponent const * unsafeGetRepresentation() const { return mValue; }


  exponent component() const { return *mValue; }

#ifdef MATHICGB_USE_MONOID
  operator Monoid::ConstMonoRef() const {
    return Monoid::toRef(mValue);
  }
#endif

private:
  const exponent& operator[](size_t index) const { return mValue[index]; }
  const exponent& operator*() const { return *mValue; }

protected:
  exponent const* mValue;
};

class Monomial : public ConstMonomial
{
  friend class PolyRing;
  friend class OrderA;
  friend class OrderB;
  friend class OrderC;
  friend class OrderD;
  friend class OrderE;
public:
  //* Wrap a raw pointer to create a monomial
  // Assumptions:
  //  1. This is done in the presence of a PolyRing
  //  2. Space for the monomial has been created
  Monomial() : ConstMonomial() {}
  Monomial(exponent *val) : ConstMonomial(val) {}

  void swap(Monomial& monomial) {
    std::swap(mValue, monomial.mValue);
  }

  exponent * unsafeGetRepresentation() { return const_cast<exponent *>(mValue); }
  exponent const * unsafeGetRepresentation() const { return mValue; }

#ifdef MATHICGB_USE_MONOID
  operator Monoid::MonoRef() {
    return Monoid::toRef(unsafeGetRepresentation());
  }
#endif

private:
  const exponent& operator[](size_t index) const { return mValue[index]; }
  exponent& operator[](size_t index) { return unsafeGetRepresentation()[index]; }

  const exponent& operator*() const { return *mValue; }
  exponent& operator*() { return * const_cast<exponent *>(mValue); }
};

inline const Monomial& ConstMonomial::castAwayConst() const
{
  return reinterpret_cast<const Monomial&>(*this);
}

typedef Monomial monomial;
typedef ConstMonomial const_monomial;

struct const_term {
  const_coefficient coeff;
  const_monomial monom; // includes component

  const_term() {}
  const_term(const_coefficient c, const_monomial m) : coeff(c), monom(m) {}
};

struct term {
  coefficient coeff;
  monomial monom; // includes component
  operator const_term() const {return const_term(coeff, monom);}

  term() {}
  term(coefficient c, monomial m) : coeff(c), monom(m) {}
};

class PolyRing {
public:
  PolyRing(coefficient charac, int nvars, const std::vector<exponent>& weights);
  PolyRing(coefficient charac, int nvars, int nweights);
  ~PolyRing() {}

  size_t getMemoryUse() const {
    // todo: Make this more accurate.
    return mMonomialPool.getMemoryUse();
  }

  coefficient charac() const { return mCharac; }
  size_t getNumVars() const { return mNumVars; }
  //       const std::vector<int> &degs,
  //       const std::string &monorder);

  static PolyRing *read(std::istream &i);
  void write(std::ostream &o) const;
  // Format for ring
  //   <char> <mNumVars> <deg1> ... <deg_n> <monorder>

  void printRingFrobbyM2Format(std::ostream& out) const;

  //  Allocate a monomial from an arena A
  //  This monomial may only be freed if no other elements that were allocated
  // later are live on A.  In this case, use freeMonomial(A,m) to free 'm'.
  Monomial allocMonomial(memt::Arena &A) const {
    exponent* ptr = static_cast<exponent*>(A.alloc(mMaxMonomialByteSize));
#ifdef MATHICGB_DEBUG
    // fill with value that do not make sense to catch bugs in debug
    // mode. The maximum value of setting all bits increases the
    // chances of getting an assert.
    std::fill_n(reinterpret_cast<char*>(ptr), mMaxMonomialByteSize,
                ~static_cast<char>(0));
#endif
    return ptr;
  }

  // Free monomial 'm' obtained by allocMonomial(A) by calling
  // freeMonomial(A,m) Recall this only works if this was the last
  // thing allocated in A.
  void freeTopMonomial(memt::Arena &A, Monomial m) const {
    A.freeTop(m.unsafeGetRepresentation());
  }

  bool fromPool(ConstMonomial m) const {
    return mMonomialPool.fromPool(m.unsafeGetRepresentation());
  }

  //  Allocate a monomial from a pool that has had its size set to 
  //   maxMonomialByteSize()
  //  Free monomials here using the SAME pool
  Monomial allocMonomial(memt::BufferPool &P) const {
    exponent* ptr = static_cast<exponent*>(P.alloc());
#ifdef MATHICGB_DEBUG
    // fill with value that do not make sense to catch bugs in debug
    // mode. The maximum value of setting all bits increases the
    // chances of getting an assert.
    std::fill_n(reinterpret_cast<char*>(ptr), mMaxMonomialByteSize,
                ~static_cast<char>(0));
#endif
    return ptr;
  }

  // Free monomial 'm' obtained by allocMonomial(P) 
  // by calling freeMonomial(P,m)
  void freeMonomial(memt::BufferPool &P, Monomial m) const {
    P.free(m.unsafeGetRepresentation());
  }

  // Only call this method for monomials returned by allocMonomial().
  void freeMonomial(Monomial m) const {
    mMonomialPool.free(m.unsafeGetRepresentation());
  }

  // Free monomials allocated here by calling freeMonomial().
  monomial allocMonomial() const { return allocMonomial(mMonomialPool); }



  coefficient toCoefficient(int64 value) const;
  coefficient coefficientNegate(coefficient coeff) const;
  coefficient coefficientNegateNonZero(coefficient coeff) const;
  coefficient coefficientSubtract(coefficient a, coefficient b) const;

  void coefficientFromInt(coefficient &result, int a) const;
  void coefficientSetOne(coefficient &result) const { result = 1; }
  void coefficientAddOneTo(coefficient &result) const;
  void coefficientReciprocalTo(coefficient &result) const;
  void coefficientNegateTo(coefficient &result) const; // result = -result
  void coefficientAddTo(coefficient &result, coefficient a, coefficient b) const; // result += a*b
  void coefficientAddTo(coefficient &result, coefficient a) const; // result += a
  void coefficientMultTo(coefficient &result, coefficient a) const; // result *= a
  void coefficientMult(coefficient a, coefficient b, coefficient &result) const; // result = a*b
  void coefficientDivide(coefficient a, coefficient b, coefficient &result) const; // result /= a
  void coefficientSet(coefficient &result, coefficient a) const { result = a; }
  bool coefficientIsZero(coefficient a) const { return a == 0; }
  bool coefficientIsOne(coefficient a) const { return a == 1; }
  bool coefficientEQ(coefficient a, coefficient b) const { return a == b; }

  // Format for monomials might be changeable, but for now, let us just do
  // array of ints:
  //  -comp exp0 exp1 ... exp(mNumVars-1) -wtr ... -wt1
  // However when written and read from a file the format will be:
  //  nterms comp v1 e1 ... v_nterms e_nterms
  // with each e_i nonzero, and v_1 > v_2 > ... > v_nterms

  size_t maxMonomialByteSize() const { return mMaxMonomialByteSize; }

  size_t maxMonomialSize() const { return mMaxMonomialSize; }

  void displayHashValues() const;

  size_t monomialHashIndex() const { return mHashIndex; }

  ///////////////////////////////////////////
  // Monomial Routines //////////////////////
  ///////////////////////////////////////////

  HashValue monomialHashValue(ConstMonomial m) const {
#ifdef MATHICGB_USE_MONOID
    return monoid().hash(m);
#else
    return static_cast<HashValue>(m[mHashIndex]);
#endif
  }

  void monomialSetExponent(Monomial m, size_t var, exponent c) const {
#ifdef MATHICGB_USE_MONOID
    monoid().setExponent(var, c, m);
#else
    m[var+1] = c;
#endif
  }

  void monomialSetExponents(Monomial m, exponent* exponents) const {
#ifdef MATHICGB_USE_MONOID
    monoid().setExponents(exponents, m);
#else
    *m = 0;
    std::memcpy(
      m.unsafeGetRepresentation() + 1,
      exponents,
      mNumVars * sizeof(exponent)
    );
    setWeightsAndHash(m);
#endif
  }

  exponent monomialExponent(ConstMonomial m, size_t var) const {
#ifdef MATHICGB_USE_MONOID
    return monoid().exponent(m, var);
#else
    return m[var+1];
#endif
  }

  // This function only sets component and the monomial itself. NOT weights, degree, or hash value
  //TODO: get Bjarke to name this function!!
  void mysteriousSPairMonomialRoutine(ConstMonomial newSig,
                                      ConstMonomial newLead,
                                      ConstMonomial baseDivSig,
                                      ConstMonomial baseDivLead,
                                      Monomial result) const;

  // Returns the weight (degree) of a. Takes the first weight if
  // working with several weight vectors.
  exponent weight(ConstMonomial a) const;

  void setWeightsAndHash(Monomial& a) const;

  inline void setWeightsOnly(Monomial& a) const;

  inline void setHashOnly(Monomial& a) const;

  bool hashValid(const_monomial m) const;

  bool weightsCorrect(ConstMonomial a) const;

  // returns LT, EQ, or GT, depending on sig ? (m2 * sig2).
  int monomialCompare(ConstMonomial a, 
                      ConstMonomial b) const; 
  // returns LT, EQ or GT
  int monomialCompare(ConstMonomial sig, 
                      ConstMonomial m2, 
                      ConstMonomial sig2) const;

  // If this method returns true for monomials a and b then it is guaranteed
  // the multiplying a and b together will not overflow the underlying
  // exponent integer. Does not work for negative exponents.
  bool monomialHasAmpleCapacity(ConstMonomial mono) const;

  bool monomialLT(ConstMonomial a, ConstMonomial b) const {
#ifdef MATHICGB_USE_MONOID
    return monoid().lessThan(a, b);
#else
    for (size_t i = mTopIndex; i != static_cast<size_t>(-1); --i)
      {
        exponent cmp = a[i] - b[i];
        if (cmp == 0) continue;
        if (cmp < 0) return false;
        return true;
      }
    return false;
#endif
  }

  bool monomialEQ(ConstMonomial a, ConstMonomial b) const;

  /// as monomialEQ, but optimized for the case that the answer is true.
  bool monomialEqualHintTrue(ConstMonomial a, ConstMonomial b) const;

  size_t monomialSize(ConstMonomial) const { return mMaxMonomialSize; }

  exponent monomialGetComponent(ConstMonomial a) const { return *a.mValue; }

  void monomialChangeComponent(Monomial a, int x) const {
#ifdef MATHICGB_USE_MONOID
    monoid().setComponent(x, a);
#else
    a[mHashIndex] -= static_cast<HashValue>(*a.mValue);
    a[mHashIndex] += static_cast<HashValue>(x);
    *a = x;
#endif
  }

  void monomialSetIdentity(Monomial& result) const;

  void monomialEi(size_t i, Monomial &result) const;

  void monomialMult(ConstMonomial a, ConstMonomial b, Monomial &result) const;

  void monomialMultTo(Monomial &a, ConstMonomial b) const; // a *= b

  /// Result is set to b/a. a must divide b.
  void monomialDivide(ConstMonomial a, ConstMonomial b, Monomial &result) const;

  /// sets result to a/b, even if that produces negative exponents.
  void monomialDivideToNegative(ConstMonomial a, ConstMonomial b, Monomial &result) const;

  /// Sets aColonB = a:b and bColonA = b:a.
  void monomialColons(ConstMonomial a, ConstMonomial b, monomial aColonB, monomial bColonA) const;

  /// returns true if b divides a.  Components are ignored.
  bool monomialIsDivisibleBy(ConstMonomial a, ConstMonomial b) const;

  /// Returns true if ab is the product of a and b.
  bool monomialIsProductOf
    (ConstMonomial a, ConstMonomial b, ConstMonomial ab) const;

  /// As monomialIsProductOf but optimized for the case that the result
  /// is true.
  bool monomialIsProductOfHintTrue
    (ConstMonomial a, ConstMonomial b, ConstMonomial ab) const;

  /// As monomialIsProductOfHintTwo(), but checks two products are equal.
  /// The return value is true if a1*b = a1b and a2*b = a2b.
  MATHICGB_INLINE bool monomialIsTwoProductsOfHintTrue(
    ConstMonomial a1,
    ConstMonomial a2,
    ConstMonomial b,
    ConstMonomial a1b,
    ConstMonomial a2b) const;

  /// Returns the hash of the product of a and b.
  HashValue monomialHashOfProduct(ConstMonomial a, ConstMonomial b) const {
    return static_cast<exponent>(
      static_cast<HashValue>(a[mHashIndex]) +
      static_cast<HashValue>(b[mHashIndex]));
  }

  void monomialCopy(ConstMonomial  a, Monomial &result) const;

  void monomialQuotientAndMult(ConstMonomial a, 
                               ConstMonomial b, 
                               ConstMonomial c, 
                               Monomial& result) const;
  // result is set to (a//b) * c

  inline bool monomialRelativelyPrime(ConstMonomial a, 
                                      ConstMonomial b) const;

  void monomialFindSignature(ConstMonomial v1,
                             ConstMonomial v2,
                             ConstMonomial u1,
                             Monomial& t1) const; 
  // answer into the already allocated t1

  size_t monomialSizeOfSupport(ConstMonomial m) const;

  void monomialGreatestCommonDivisor(ConstMonomial a, 
                                     ConstMonomial b, 
                                     Monomial& g) const;

  inline void monomialLeastCommonMultiple(ConstMonomial a, 
                                          ConstMonomial b, 
                                          Monomial& l) const;

  inline void monomialLeastCommonMultipleNoWeights(ConstMonomial a, 
                                                   ConstMonomial b, 
                                                   Monomial& l) const;

  bool monomialIsLeastCommonMultiple(ConstMonomial a, 
                                     ConstMonomial b, 
                                     ConstMonomial l) const;

  bool monomialIsLeastCommonMultipleNoWeights(ConstMonomial a, 
                                              ConstMonomial b, 
                                              ConstMonomial l) const;

  // Returns true if there is a variable var such that hasLarger raises var to
  // a strictly greater exponent than both smaller1 and smaller2 does.
  inline bool monomialHasStrictlyLargerExponent(
    ConstMonomial hasLarger,
    ConstMonomial smaller1,
    ConstMonomial smaller2) const;

  void monomialParse(std::istream& i, 
                     Monomial& result) const;

  void monomialDisplay(std::ostream& o, 
                       ConstMonomial a, 
                       bool print_comp=true, 
                       bool print_one=true) const;
  void monomialDisplay(FILE* file,
                       ConstMonomial a, 
                       bool printComponent = true, 
                       bool printOne = true) const;

  void printMonomialFrobbyM2Format(std::ostream& out, ConstMonomial m) const;

  ///////////////////////////////////////////
  ///////////////////////////////////////////

  struct coefficientStats {
    size_t n_addmult;
    size_t n_add;
    size_t n_mult;
    size_t n_recip;
    size_t n_divide;
  };
  const coefficientStats & getCoefficientStats() const { return mStats; }
  void resetCoefficientStats() const;

private:
  inline HashValue computeHashValue(const_monomial a1) const;

  coefficient mCharac; // p=mCharac: ring is ZZ/p
  size_t mNumVars;
  size_t mNumWeights; // stored as negative of weight vectors
  size_t mTopIndex;
  size_t mHashIndex; // 1 more than mTopIndex.  Where the has value is stored.
  size_t mMaxMonomialSize;
  size_t mMaxMonomialByteSize;
  std::vector<exponent> mWeights; // 0..mNumWeights * mNumVars - 1.

  std::vector<HashValue> mHashVals; // one for each variable 0..mNumVars-1
  // stored as weightvec1 weightvec2 ...

  mutable memt::BufferPool mMonomialPool;
  mutable coefficientStats mStats;

  bool mTotalDegreeGradedOnly;

#ifdef MATHICGB_USE_MONOID
  const Monoid& monoid() const {return mMonoid;}
  Monoid mMonoid;
#endif
#ifdef MATHICGB_USE_FIELD
  const Field field() const {return mField;}
  Field mField;
#endif
};

inline exponent PolyRing::weight(ConstMonomial a) const {
#ifdef MATHICGB_USE_MONOID
  return monoid().degree(a);
#else
  MATHICGB_ASSERT(weightsCorrect(a));
  return a[mNumVars + 1];
#endif
}

////////////////////////////////////////////////
// New Monomial Routines ///////////////////////
////////////////////////////////////////////////

inline bool PolyRing::monomialEQ(ConstMonomial a, ConstMonomial b) const
{
#ifdef MATHICGB_USE_MONOID
  return monoid().equal(a, b);
#else
  for (size_t i = 0; i <= mNumVars; ++i)
    if (a[i] != b[i])
      return false;
  return true;
#endif
}

inline bool PolyRing::monomialEqualHintTrue(
  const ConstMonomial a,
  const ConstMonomial b
) const {
#ifdef MATHICGB_USE_MONOID
  return monoid().equalHintTrue(a, b);
#else
  // if a[i] != b[i] then a[i] ^ b[i] != 0, so the or of all xors is zero
  // if and only if a equals b. This way we avoid having a branch to check
  // equality for every iteration of the loop, which is a win in the case
  // that none of the early-exit branches are taken - that is, when a equals b.
  exponent orOfXor = 0;
  for (size_t i = mNumVars; i != 0; --i)
    orOfXor |= a[i] ^ b[i];
  const bool areEqual = (orOfXor == 0);
  MATHICGB_ASSERT(areEqual == monomialEQ(a, b));
  return areEqual;
#endif
}

inline bool PolyRing::monomialIsProductOfHintTrue(
  const ConstMonomial a, 
  const ConstMonomial b, 
  const ConstMonomial ab
) const {
#ifdef MATHICGB_USE_MONOID
  return monoid().isProductOfHintTrue(a, b, ab);
#else
  // We compare more than one exponent at a time using 64 bit integers. This 
  // might go one 32 bit value at the end too far, but since that space is
  // either a degree or a hash value that is fine --- those values will also
  // match if the monomials are equal. This does not work for negative
  // exponents since the overflowing bit will go into the next word.
  // It is OK that the degree field can be negative (a field we might go
  // into without caring about it because it shares a 64 bit field with
  // the last exponent), because it is at the end so the overflowing
  // bit will not interfere.

  // todo: ensure 8 byte alignment. Though there seem to be no ill effects
  // for unaligned access. Performance seems to be no worse than for using
  // 32 bit integers directly.

  if (sizeof(exponent) < 4)
    return monomialIsProductOf(a, b, ab);

  uint64 orOfXor = 0;
  for (size_t i = mNumVars / 2; i != static_cast<size_t>(-1); --i) {
    MATHICGB_ASSERT(a[i*2] >= 0);
    MATHICGB_ASSERT(i == mNumVars / 2 || a[i*2+1] >= 0);

    uint64 A, B, AB;
    // We have to use std::memcpy here because just casting to a int64 breaks
    // the strict aliasing rule which implies undefined behavior. Both MSVC and
    // gcc don't actually call memcpy here. MSVC is a tiny bit slower for this
    // code than for casting while GCC seems to be exactly the same speed.
    std::memcpy(&A, &a[i*2], 8);
    std::memcpy(&B, &b[i*2], 8);
    std::memcpy(&AB, &ab[i*2], 8);
    orOfXor |= AB ^ (A + B);
  }
  MATHICGB_ASSERT((orOfXor == 0) == monomialIsProductOf(a, b, ab));

  return orOfXor == 0; 
#endif
}

MATHICGB_INLINE bool PolyRing::monomialIsTwoProductsOfHintTrue(
  const ConstMonomial a1,
  const ConstMonomial a2,
  const ConstMonomial b,
  const ConstMonomial a1b,
  const ConstMonomial a2b
) const {
#ifdef MATHICGB_USE_MONOID
  return monoid().isTwoProductsOfHintTrue(a1, a2, b, a1b, a2b);
#else
  if (sizeof(exponent) < 4)
    return (monomialIsProductOf(a1, b, a1b) &&
      monomialIsProductOf(a2, b, a2b));

  uint64 orOfXor = 0;
  for (size_t i = mNumVars / 2; i != static_cast<size_t>(-1); --i) {
    uint64 A1, A2, B, A1B, A2B;
    std::memcpy(&A1, &a1[i*2], 8);
    std::memcpy(&A2, &a2[i*2], 8);
    std::memcpy(&B, &b[i*2], 8);
    std::memcpy(&A1B, &a1b[i*2], 8);
    std::memcpy(&A2B, &a2b[i*2], 8);
    orOfXor |= (A1B ^ (A1 + B)) | (A2B ^ (A2 + B));
  }
  MATHICGB_ASSERT((orOfXor == 0) ==
    (monomialIsProductOf(a1, b, a1b) && monomialIsProductOf(a2, b, a2b)));

  return orOfXor == 0;
#endif
}

inline bool PolyRing::monomialIsProductOf(
  ConstMonomial a, 
  ConstMonomial b, 
  ConstMonomial ab
) const {
#ifdef MATHICGB_USE_MONOID
  return monoid().isProductOf(a, b, ab);
#else
  for (size_t i = 0; i <= mNumVars; ++i)
    if (ab[i] != a[i] + b[i])
      return false;
  return true;
#endif
}

inline void PolyRing::monomialMult(ConstMonomial a, 
                                   ConstMonomial b, 
                                   Monomial &result) const
{
#ifdef MATHICGB_USE_MONOID
  monoid().multiply(a, b, result);
#else
  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    result[i] = a[i] + b[i];
  MATHICGB_ASSERT(computeHashValue(result) ==
    static_cast<exponent>(computeHashValue(a) + computeHashValue(b)));

#if 0
  // testing different things to see if we can speed it up further.
  // changing to ascending loop slowed it down.
  // ascending, with pointers: slightly faster than prev, but still slower than above simple code
  exponent *presult = result.unsafeGetRepresentation();
  exponent const * pa = a.unsafeGetRepresentation();
  exponent const * pb = b.unsafeGetRepresentation();
  for (size_t i=0; i<= mHashIndex; ++i)
    //  for (size_t i = mHashIndex; i != static_cast<size_t>(-1); --i)
    *presult++ = *pa++ + *pb++;
  //    result[i] = a[i] + b[i];
#endif
#endif
}

inline void PolyRing::setWeightsOnly(Monomial& a1) const
{
  exponent *a = a1.unsafeGetRepresentation();
  a++;
  auto wts = mWeights.data();
  for (size_t i = 0; i < mNumWeights; ++i)
    {
      exponent result = 0;
      for (size_t j = 0; j < mNumVars; ++j)
        result += *wts++ * a[j];
      a[mNumVars+i] = result;
    }
}

inline HashValue PolyRing::computeHashValue(const_monomial a1) const {
#ifdef MATHICGB_USE_MONOID
  return monoid().computeHash(a1);
#else
  const exponent* a = a1.unsafeGetRepresentation();
  HashValue hash = static_cast<HashValue>(*a);
  a++;
  for (size_t i = 0; i < mNumVars; ++i)
    hash += static_cast<HashValue>(a[i]) * mHashVals[i];
  // cast to potentially discard precision that will also be lost
  // when storing a hash value as an exponent. Otherwise the hash
  // value that is computed will not match the stored hash value.
  return static_cast<exponent>(hash);
#endif
}

inline void PolyRing::setHashOnly(Monomial& a1) const
{
  exponent* a = a1.unsafeGetRepresentation();
  a[mHashIndex] = computeHashValue(a1);
}

inline int PolyRing::monomialCompare(ConstMonomial a, ConstMonomial b) const
// returns LT, EQ or GT
{
#ifdef MATHICGB_USE_MONOID
  return monoid().compare(a, b);
#else
  for (size_t i = mTopIndex; i != static_cast<size_t>(-1); --i)
    {
      auto cmp = a[i] - b[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
  return EQ;
#endif
}

inline bool PolyRing::monomialIsDivisibleBy(ConstMonomial a,
                                            ConstMonomial b) const
{
#ifdef MATHICGB_USE_MONOID
  return monoid().divides(b, a);
#else
  for (size_t i = 1; i<= mNumVars; i++)
    if (a[i] < b[i])
      return false;
  return true;
#endif
}

inline void PolyRing::monomialDivide(ConstMonomial a, 
                                     ConstMonomial b, 
                                     Monomial& result) const
{
#ifdef MATHICGB_USE_MONOID
  return monoid().divide(b, a, result);
#else
  //// returns true if b divides a, in this case, result is set to b//a.
  size_t i;
  for (i = 1; i <= mNumVars; i++)
    {
      exponent c = a[i] - b[i];
      if (c < 0)
        return;
      result[i] = c;
    }
  // at this point we have divisibility, so need to fill in the rest of the monomial
  *result = *a.mValue - *b.mValue;  // component
  for ( ; i<=mHashIndex; i++)
    result[i] = a[i] - b[i];
#endif
}

inline void PolyRing::monomialColons(
  ConstMonomial a,
  ConstMonomial b,
  monomial aColonB,
  monomial bColonA
) const {
#ifdef MATHICGB_USE_MONOID
  monoid().colons(a, b, aColonB, bColonA);
#else
  *aColonB = *a;
  *bColonA = *b;
  for (size_t i = 1; i <= mNumVars; i++) {
    exponent max = std::max(a[i], b[i]);
    aColonB[i] = max - b[i];
    bColonA[i] = max - a[i];
  }
  setWeightsAndHash(aColonB);
  setWeightsAndHash(bColonA);
#endif
}

inline void PolyRing::monomialDivideToNegative(ConstMonomial a, 
                                               ConstMonomial b, 
                                               Monomial& result) const 
{
#ifdef MATHICGB_USE_MONOID
  monoid().divideToNegative(b, a, result);
#else
  for (size_t i = 0; i <= mHashIndex; ++i)
    result[i] = a[i] - b[i];
  MATHICGB_ASSERT(monomialHashValue(result) ==
    static_cast<exponent>(monomialHashValue(a) - monomialHashValue(b)));
  MATHICGB_ASSERT(!hashValid(a) || !hashValid(b) || hashValid(result));
  MATHICGB_ASSERT(computeHashValue(result) == static_cast<exponent>
    (computeHashValue(a) - computeHashValue(b)));
#endif
}

inline bool PolyRing::monomialRelativelyPrime(ConstMonomial a, 
                                              ConstMonomial b) const
{
#ifdef MATHICGB_USE_MONOID
  return monoid().relativelyPrime(a, b);
#else
  for (size_t i = 1; i <= mNumVars; ++i)
    if (a[i] > 0 && b[i] > 0)
      return false;
  return true;
#endif
}

inline void PolyRing::monomialLeastCommonMultiple(
  ConstMonomial a,
  ConstMonomial b,
  Monomial& l) const
{
#ifdef MATHICGB_USE_MONOID
  monoid().lcm(a, b, l);
#else
  monomialLeastCommonMultipleNoWeights(a, b, l);
  setWeightsAndHash(l);
#endif
}

inline void PolyRing::monomialLeastCommonMultipleNoWeights(
  ConstMonomial a,
  ConstMonomial b,
  Monomial& l) const
{
#ifdef MATHICGB_USE_MONOID
  monoid().lcmRaw(a, b, l);
#else
  *l = 0;
  for (size_t i = 1; i <= mNumVars; ++i)
    l[i] = std::max(a[i], b[i]);
#endif
}

inline bool PolyRing::monomialHasStrictlyLargerExponent(
  ConstMonomial hasLarger,
  ConstMonomial smaller1,
  ConstMonomial smaller2) const 
{
#ifdef MATHICGB_USE_MONOID
  return !monoid().dividesLcm(hasLarger, smaller1, smaller2);
#else
  for (size_t i = 1; i <= mNumVars; ++i)
    if (hasLarger[i] > smaller1[i] && hasLarger[i] > smaller2[i])
      return true;
  return false;
#endif
}


////////////////////////////////////////////////
// Old Monomial Routines ///////////////////////
////////////////////////////////////////////////

inline bool PolyRing::monomialIsLeastCommonMultiple(
  ConstMonomial a,
  ConstMonomial b,
  ConstMonomial l) const
{
#ifdef MATHICGB_USE_MONOID
  return monoid().isLcm(a, b, l);
#else
  return monomialIsLeastCommonMultipleNoWeights(a, b, l) && weightsCorrect(l);
#endif
}

inline void PolyRing::coefficientReciprocalTo(coefficient& result) const
{
#ifdef MATHICGB_USE_FIELD
  result = field().inverse(field().toElementInRange(result)).value();
#else
  MATHICGB_ASSERT(result != 0);
  mStats.n_recip++;
  result = modularInverse(result, mCharac);
#endif
}

inline void PolyRing::coefficientDivide(coefficient a, coefficient b, coefficient &result) const
 // result = a/b
{
#ifdef MATHICGB_USE_FIELD
  result = field().quotient
    (field().toElementInRange(a), field().toElementInRange(b)).value();
#else
  mStats.n_divide++;
  result = (a * modularInverse(b, mCharac)) % mCharac;
  MATHICGB_ASSERT((result * b) % mCharac == a);
  MATHICGB_ASSERT(result >= 0);
  MATHICGB_ASSERT(result < mCharac);
#endif
}

inline void PolyRing::coefficientFromInt(coefficient &result, int a) const
{
#ifdef MATHICGB_USE_FIELD
  result = field().toElement(a).value();
#else
  result = toCoefficient(a);
#endif
}

inline void PolyRing::coefficientAddOneTo(coefficient &result) const
{
#ifdef MATHICGB_USE_FIELD
  result = field().plusOne(field().toElementInRange(result)).value();
#else
  ++result;
  if (result == mCharac)
    result = 0;
#endif
}

inline void PolyRing::coefficientNegateTo(coefficient& result) const {
#ifdef MATHICGB_USE_FIELD
  result = field().negative(field().toElementInRange(result)).value();
#else
  MATHICGB_ASSERT(result < mCharac);
  if (result != 0)
    result = coefficientNegateNonZero(result);
#endif
}

inline coefficient PolyRing::toCoefficient(const int64 value) const {
#ifdef MATHICGB_USE_FIELD
  return field().toElement(value).value();
#else
  auto modLong = value % mCharac;
  if (modLong < 0)
    modLong += mCharac;
  MATHICGB_ASSERT(0 <= modLong);
  MATHICGB_ASSERT(modLong < mCharac);
  const auto mod = static_cast<coefficient>(modLong);
  MATHICGB_ASSERT(0 <= mod);
  MATHICGB_ASSERT(mod < mCharac);
  return mod;
#endif
}

inline coefficient PolyRing::coefficientNegate(const coefficient coeff) const {
#ifdef MATHICGB_USE_FIELD
  return field().negative(field().toElementInRange(coeff)).value();
#else
  MATHICGB_ASSERT(coeff < mCharac);
  return coeff == 0 ? 0 : coefficientNegateNonZero(coeff);
#endif
}

inline coefficient PolyRing::coefficientNegateNonZero(
  const coefficient coeff
) const {
#ifdef MATHICGB_USE_FIELD
  return field().negativeNonZero(field().toElementInRange(coeff)).value();
#else
  MATHICGB_ASSERT(coeff != 0);
  MATHICGB_ASSERT(coeff < mCharac);
  return mCharac - coeff;
#endif
}

inline coefficient PolyRing::coefficientSubtract(
  const coefficient a,
  const coefficient b
) const {
#ifdef MATHICGB_USE_FIELD
  return field().difference
    (field().toElementInRange(a), field().toElementInRange(b)).value();
#else
  MATHICGB_ASSERT(a < mCharac);
  MATHICGB_ASSERT(b < mCharac);
  const auto diff = a < b ? a + (mCharac - b) : a - b;
  MATHICGB_ASSERT(diff < mCharac);
  MATHICGB_ASSERT((diff + b) % mCharac == a);
  return diff;
#endif
}

inline void PolyRing::coefficientAddTo
(coefficient &result, coefficient a, coefficient b) const
// result += a*b
{
#ifdef MATHICGB_USE_FIELD
  const auto prod =
    field().product(field().toElementInRange(a), field().toElementInRange(b));
  result = field().sum(field().toElementInRange(result), prod).value();
#else
  mStats.n_addmult++;
  auto c = a * b + result;
  result = c % mCharac;
#endif
}

inline void PolyRing::coefficientAddTo(coefficient &result, coefficient a) const
 // result += a
{
#ifdef MATHICGB_USE_FIELD
  result = field().sum
    (field().toElementInRange(result), field().toElementInRange(a)).value();
#else
  mStats.n_add++;
  result += a;
  if (result >= mCharac)
    result -= mCharac;
#endif
}

inline void PolyRing::coefficientMultTo
(coefficient &result, coefficient a) const
  // result *= a
{
#ifdef MATHICGB_USE_FIELD
  result = field().product
    (field().toElementInRange(result), field().toElementInRange(a)).value();
#else
  mStats.n_mult++;
  coefficient b = result * a;
  result = b % mCharac;
#endif
}

inline void PolyRing::coefficientMult
(coefficient a, coefficient b, coefficient &result) const
{
#ifdef MATHICGB_USE_FIELD
  result = field().product
    (field().toElementInRange(a), field().toElementInRange(b)).value();
#else
  mStats.n_mult++;
  coefficient c = b * a;
  result = c % mCharac;
#endif
}

inline bool PolyRing::monomialHasAmpleCapacity(ConstMonomial mono) const {
#ifdef MATHICGB_USE_MONOID
  return monoid().hasAmpleCapacity(mono);
#else
  const auto halfMax = std::numeric_limits<exponent>::max() / 2;
  for (size_t i = mTopIndex; i != 0; --i)
    if (mono[i] > halfMax)
      return false;
  return true;
#endif
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
