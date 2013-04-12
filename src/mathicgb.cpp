#include "mathicgb/stdinc.h"
#include "mathicgb.h"

#include "mathicgb/Ideal.hpp"
#include "mathicgb/PolyRing.hpp"
#include "mathicgb/Poly.hpp"
#include "mathicgb/Reducer.hpp"
#include "mathicgb/BuchbergerAlg.hpp"
#include <mathic.h>

namespace {
  bool isPrime(unsigned int n) {
    if (n == 0 || n == 1)
      return false;
    if (n == 2 || n == 3)
      return true;
    return true; // todo: make better test
  }
};

#ifndef MATHICGB_ASSERT
#ifdef MATHICGB_DEBUG
#include <cassert>
#define MATHICGB_ASSERT(X) assert(X)
#else
#define MATHICGB_ASSERT(X)
#endif
#endif

#define MATHICGB_STREAM_CHECK(X, MSG) \
  do { \
    const bool value = (X); \
    if (!value) { \
      const bool ignoreMe = false; \
      MATHICGB_ASSERT(( \
        "MathicGB stream protocol error: "#MSG \
        "\nAssert expression: "#X"\n", \
        false \
      )); \
      throw std::invalid_argument( \
        "MathicGB stream protocol error: "#MSG \
        "\nAssert expression: "#X"\n" \
      ); \
    } \
  } while (false)

namespace mgbi {
  struct StreamStateChecker::Pimpl {
    Pimpl(Coefficient modulus, VarIndex varCount):
      modulus(modulus),
      varCount(varCount),
      state(Initial),

      hasClaimedPolyCount(false),
      claimedPolyCount(0),
      seenPolyCount(0),

      hasClaimedTermCount(false),
      claimedTermCount(0),
      seenTermCount(0),

      lastVar(0)
    {}

    bool debugAssertValid() const;

    enum State {
      Initial,
      MakingIdeal,
      MakingPoly,
      MakingTerm,
      HasIdeal
    };

    const Coefficient modulus;
    const VarIndex varCount;

    State state;

    bool hasClaimedPolyCount;
    size_t claimedPolyCount;
    size_t seenPolyCount;

    bool hasClaimedTermCount;
    size_t claimedTermCount;
    size_t seenTermCount;

    VarIndex lastVar;
  };

  bool StreamStateChecker::Pimpl::debugAssertValid() const {
#ifdef MATHICGB_DEBUG
    MATHICGB_ASSERT(this != 0);
    switch (state) {
    case Initial:
    case MakingIdeal:
    case MakingPoly:
    case MakingTerm:
    case HasIdeal:
      break;

    default:
      MATHICGB_ASSERT(false);
      return false;
    }
#endif
    return true;
  };

  StreamStateChecker::StreamStateChecker(const Coefficient modulus, const VarIndex varCount):
    mPimpl(new Pimpl(modulus, varCount))
  {
    try {
      MATHICGB_STREAM_CHECK(isPrime(modulus), "The modulus must be prime");
      MATHICGB_ASSERT(mPimpl->debugAssertValid());
    } catch (...) {
      delete mPimpl;
    }
  }

  StreamStateChecker::~StreamStateChecker() {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());
    delete mPimpl;
  }

  void StreamStateChecker::idealBegin() {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());

    MATHICGB_STREAM_CHECK(
      mPimpl->state == Pimpl::Initial || mPimpl->state == Pimpl::HasIdeal,
      "idealBegin() must not be called twice "
      "without an intervening call to idealDone()."
    );
    mPimpl->state = Pimpl::MakingIdeal;
    mPimpl->hasClaimedPolyCount = false;
    mPimpl->claimedPolyCount = 0;
    mPimpl->seenPolyCount = 0;

    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::idealBegin(size_t polyCount) {
    idealBegin();
    mPimpl->hasClaimedPolyCount = true;
    mPimpl->claimedPolyCount = polyCount;
    mPimpl->seenPolyCount = 0;

    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::appendPolynomialBegin() {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());

    MATHICGB_STREAM_CHECK(
      mPimpl->state != Pimpl::Initial && mPimpl->state != Pimpl::HasIdeal,
      "appendPolynomialBegin() must only be called after idealBegin() and "
      "before idealEnd()."
    );
    MATHICGB_STREAM_CHECK(
      mPimpl->state == Pimpl::MakingIdeal,
      "appendPolynomialBegin() must not be called twice without "
      "an intervening call to appendPolynomialDone()."
    );
    MATHICGB_STREAM_CHECK(
      !mPimpl->hasClaimedPolyCount ||
        mPimpl->seenPolyCount < mPimpl->claimedPolyCount,
      "The number of polynomials in an ideal must not exceed the amount "
      "passed to idealBegin()."
    );
    mPimpl->state = Pimpl::MakingPoly;
    mPimpl->seenPolyCount += 1;
    mPimpl->hasClaimedTermCount = false;
    mPimpl->claimedTermCount = 0;
    mPimpl->seenTermCount = 0;

    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::appendPolynomialBegin(size_t termCount) {
    appendPolynomialBegin();
    mPimpl->hasClaimedTermCount = true;
    mPimpl->claimedTermCount = termCount;
    mPimpl->seenTermCount = 0;

    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::appendTermBegin() {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());

    MATHICGB_STREAM_CHECK(
      mPimpl->state != Pimpl::Initial &&
        mPimpl->state != Pimpl::HasIdeal &&
        mPimpl->state != Pimpl::MakingIdeal,
      "appendTermBegin() must only be called after appendPolynomialBegin() "
      "and before appendPolynomialDone()."
    );
    MATHICGB_STREAM_CHECK(
      mPimpl->state == Pimpl::MakingPoly,
      "appendTermBegin() must not be called twice without an intervening "
      "call to appendTermDone()."
    );
    MATHICGB_STREAM_CHECK(
      !mPimpl->hasClaimedTermCount ||
        mPimpl->seenTermCount < mPimpl->claimedTermCount,
      "The number of terms in a polynomial must not exceed the amount "
      "passed to appendPolynomialBegin()."
    );
     
    mPimpl->state = Pimpl::MakingTerm;
    mPimpl->seenTermCount += 1;
    mPimpl->lastVar = std::numeric_limits<decltype(mPimpl->lastVar)>::max();

    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::appendExponent(VarIndex index, Exponent exponent) {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());

    MATHICGB_STREAM_CHECK(
      mPimpl->state == Pimpl::MakingTerm,
      "appendExponent must only be called after appendTermBegin() and before "
      "appendTermDone()."
    );
    MATHICGB_STREAM_CHECK(
      index < mPimpl->varCount,
      "The index passed to appendExponent must be strictly less than "
      "the number of variables."
    );
    MATHICGB_STREAM_CHECK(
      mPimpl->lastVar ==
        std::numeric_limits<decltype(mPimpl->lastVar)>::max() ||
      mPimpl->lastVar < index,
      "The variable indices passed to appendExponent must be in strictly "
      "increasing order within each monomial."
    );

    mPimpl->lastVar = index;
    
    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::appendTermDone(Coefficient coefficient) {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());

    MATHICGB_STREAM_CHECK(
      coefficient > 0,
      "The coefficient passed to appendTermDone() must be strictly positive."
    );
    MATHICGB_STREAM_CHECK(
      coefficient < mPimpl->modulus,
      "The coefficient passed to appendTermDone() must be strictly less "
      "then the modulus."
    );
    MATHICGB_STREAM_CHECK(
      mPimpl->state == Pimpl::MakingTerm,
      "appendTermDone() must only be called after appendTermBegin()."
    );
    mPimpl->state = Pimpl::MakingPoly;
    
    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::appendPolynomialDone() {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());

    MATHICGB_STREAM_CHECK(
      mPimpl->state == Pimpl::MakingPoly,
      "appendPolynomialDone() must only be called after appendPolynomialBegin()."
    );
    MATHICGB_STREAM_CHECK(
      !mPimpl->hasClaimedTermCount ||
        mPimpl->seenTermCount == mPimpl->claimedTermCount,
      "The number of terms in a polynomial must match the amount "
      "passed to appendPolynomialBegin()."
    );
    mPimpl->state = Pimpl::MakingIdeal;

    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  void StreamStateChecker::idealDone() {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());

    MATHICGB_STREAM_CHECK(
      mPimpl->state == Pimpl::MakingIdeal,
      "idealDone() must only be called after idealBegin()."
    );
    MATHICGB_STREAM_CHECK(
      !mPimpl->hasClaimedPolyCount ||
        mPimpl->seenPolyCount == mPimpl->claimedPolyCount,
      "The number of polynomials in an ideal must match the amount "
      "passed to idealBegin()."
    );

    mPimpl->state = Pimpl::HasIdeal;

    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  bool StreamStateChecker::hasIdeal() const {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());
    return mPimpl->state == Pimpl::HasIdeal;
  }
}

namespace mgb {
  // ** Implementation of the class GroebnerConfiguration

  struct GroebnerConfiguration::Pimpl {
    Pimpl(Coefficient modulus, VarIndex varCount):
      mModulus(modulus),
      mVarCount(varCount),
      mReducer(DefaultReducer)
#ifdef MATHICGB_DEBUG
      , mHasBeenDestroyed(false)
#endif
    {
    }

    ~Pimpl() {
      MATHICGB_ASSERT(debugAssertValid());
      MATHICGB_IF_DEBUG(mHasBeenDestroyed = true;)
    }

    static bool reducerValid(Reducer reducer) {
      return
        reducer == DefaultReducer ||
        reducer == ClassicReducer ||
        reducer == MatrixReducer;
    }

    bool debugAssertValid() const {
#ifdef MATHICGB_DEBUG
      MATHICGB_ASSERT(this != 0);
      MATHICGB_ASSERT(reducerValid(mReducer));
      MATHICGB_ASSERT(mModulus != 0);
      MATHICGB_ASSERT_NO_ASSUME(!mHasBeenDestroyed);
#endif
      return true;
    }

    const Coefficient mModulus;
    const VarIndex mVarCount;
    Reducer mReducer;
    MATHICGB_IF_DEBUG(bool mHasBeenDestroyed);
  };

  GroebnerConfiguration::GroebnerConfiguration(
    Coefficient modulus,
    VarIndex varCount
  ):
    mPimpl(new Pimpl(modulus, varCount))
  {
    if (modulus > std::numeric_limits<unsigned short>::max()) {
      MATHICGB_ASSERT_NO_ASSUME(false);
      std::ostringstream str;
      str << "Modulus " << modulus
        << " is too large. MathicGB only supports 16 bit moduli.";
      mathic::reportError(str.str());
    }
    if (!isPrime(modulus)) {
      MATHICGB_ASSERT_NO_ASSUME(false);
      std::ostringstream str;
      str << "Modulus " << modulus
        << " is not prime. MathicGB only supports prime fields.";
      mathic::reportError(str.str());
    }
    MATHICGB_ASSERT(mPimpl->debugAssertValid());
  }

  GroebnerConfiguration::GroebnerConfiguration(
    const GroebnerConfiguration& conf
  ):
    mPimpl(new Pimpl(*conf.mPimpl))
  {
    MATHICGB_ASSERT(conf.mPimpl->debugAssertValid());
  }

  GroebnerConfiguration::~GroebnerConfiguration() {
    MATHICGB_ASSERT(mPimpl->debugAssertValid());
    delete mPimpl;
  }

  auto GroebnerConfiguration::modulus() const -> Coefficient {
    return mPimpl->mModulus;
  }
 
  auto GroebnerConfiguration::varCount() const -> VarIndex {
    return mPimpl->mVarCount;
  }

  void GroebnerConfiguration::setReducer(Reducer reducer) {
    MATHICGB_ASSERT(Pimpl::reducerValid(reducer));
    mPimpl->mReducer = reducer;
  }

  auto GroebnerConfiguration::reducer() const -> Reducer {
    return mPimpl->mReducer;
  }
}

// ** Implementation of class GroebnerInputIdealStream
namespace mgb {
  struct GroebnerInputIdealStream::Pimpl {
    Pimpl(const GroebnerConfiguration& conf):
      // @todo: varCount should not be int. Fix PolyRing constructor,
      // then remove this static_cast.
      ring(conf.modulus(), static_cast<int>(conf.varCount()), 1),
      ideal(ring),
      poly(ring),
      monomial(ring.allocMonomial()),
      conf(conf)
#ifdef MATHICGB_DEBUG
      , hasBeenDestroyed(false),
      checker(conf.modulus(), conf.varCount())
#endif
    {}

    ~Pimpl() {
      ring.freeMonomial(monomial);
    }

    const PolyRing ring;
    Ideal ideal;
    Poly poly;
    Monomial monomial;
    const GroebnerConfiguration conf;
    MATHICGB_IF_DEBUG(bool hasBeenDestroyed);
    MATHICGB_IF_DEBUG(::mgbi::StreamStateChecker checker); 
 };

  GroebnerInputIdealStream::GroebnerInputIdealStream(
    const GroebnerConfiguration& conf
  ):
    mExponents(new Exponent[conf.varCount()]),
    mPimpl(new Pimpl(conf))
  {
    MATHICGB_ASSERT(debugAssertValid());
  }

  GroebnerInputIdealStream::~GroebnerInputIdealStream() {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_ASSERT(mExponents != 0);
    MATHICGB_ASSERT(mPimpl != 0);
    MATHICGB_ASSERT_NO_ASSUME(!mPimpl->hasBeenDestroyed);
    MATHICGB_IF_DEBUG(mPimpl->hasBeenDestroyed = true);
    delete mPimpl;
    delete[] mExponents;
  }

  const GroebnerConfiguration& GroebnerInputIdealStream::configuration() {
    return mPimpl->conf;
  }

  auto GroebnerInputIdealStream::modulus() const -> Coefficient {
    return mPimpl->conf.modulus();
  }

  auto GroebnerInputIdealStream::varCount() const -> VarIndex {
    return mPimpl->conf.varCount();
  }

  void GroebnerInputIdealStream::idealBegin() {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.idealBegin());
    MATHICGB_ASSERT(mPimpl->poly.isZero());
    MATHICGB_ASSERT(mPimpl->ideal.empty());

    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::idealBegin(size_t polyCount) {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.idealBegin(polyCount));
    MATHICGB_ASSERT(mPimpl->poly.isZero());
    MATHICGB_ASSERT(mPimpl->ideal.empty());

    mPimpl->ideal.reserve(polyCount);

    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::appendPolynomialBegin() {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.appendPolynomialBegin());
    MATHICGB_ASSERT(mPimpl->poly.isZero());
    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::appendPolynomialBegin(size_t termCount) {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.appendPolynomialBegin(termCount));
    MATHICGB_ASSERT(mPimpl->poly.isZero());

    mPimpl->poly.reserve(termCount);

    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::appendTermBegin() {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.appendTermBegin());

    std::fill_n(mExponents, varCount(), 0);

    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::appendTermDone(Coefficient coefficient) {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.appendTermDone(coefficient));

    // @todo: do this directly into the polynomial instead of copying a second
    // time.
    mPimpl->ring.monomialSetExponents(mPimpl->monomial, mExponents);
    mPimpl->poly.appendTerm(coefficient, mPimpl->monomial);

    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::appendPolynomialDone() {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.appendPolynomialDone());

    // todo: avoid copy here by the following changes
    // todo: give Ideal a Poly&& insert
    // todo: give Poly a Poly&& constructor
    auto poly = make_unique<Poly>(std::move(mPimpl->poly));
    if (!poly->termsAreInDescendingOrder())
      poly->sortTermsDescending();
    mPimpl->ideal.insert(std::move(poly));
    mPimpl->poly.setToZero();

    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::idealDone() {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.idealDone());
  }

  bool GroebnerInputIdealStream::debugAssertValid() const {
    MATHICGB_ASSERT(this != 0);
    MATHICGB_ASSERT(mExponents != 0);
    MATHICGB_ASSERT(mPimpl != 0);
    MATHICGB_ASSERT_NO_ASSUME(!mPimpl->hasBeenDestroyed);
    MATHICGB_ASSERT(!mPimpl->monomial.isNull());
    MATHICGB_ASSERT(&mPimpl->ideal.ring() == &mPimpl->ring);
    MATHICGB_ASSERT(&mPimpl->poly.ring() == &mPimpl->ring);
    MATHICGB_ASSERT(mPimpl->ring.getNumVars() == mPimpl->conf.varCount());
    MATHICGB_ASSERT(mPimpl->ring.charac() == mPimpl->conf.modulus());
    return true;
  }
}

// ** Implementation of class mgbi::PimplOf
namespace mgbi {
  class PimplOf {
  public:
    template<class T>
    typename T::Pimpl& operator()(T& t) {
      MATHICGB_ASSERT(t.mPimpl != 0);
      return *t.mPimpl;
    }
  };
}

// ** Implementation of mgbi::IdealAdapter
namespace mgbi {
  struct IdealAdapter::Pimpl {
    std::unique_ptr<Ideal> ideal;
  };

  IdealAdapter::IdealAdapter():
    mPimpl(new Pimpl()) {
  }

  IdealAdapter::~IdealAdapter() {
    MATHICGB_ASSERT(mPimpl != 0);
    delete mPimpl;
  }

  auto IdealAdapter::varCount() const -> VarIndex {
    MATHICGB_ASSERT(mPimpl->ideal.get() != 0);
    return mPimpl->ideal->ring().getNumVars();
  }

  size_t IdealAdapter::polyCount() const {
    MATHICGB_ASSERT(mPimpl->ideal.get() != 0);
    return mPimpl->ideal->size();
  }

  size_t IdealAdapter::termCount(PolyIndex poly) const {
    MATHICGB_ASSERT(mPimpl->ideal.get() != 0);
    MATHICGB_ASSERT(poly < mPimpl->ideal->size());
    return mPimpl->ideal->getPoly(poly)->nTerms();
  }

  auto IdealAdapter::term(
    PolyIndex poly,
    TermIndex term
  ) const -> std::pair<Coefficient, const Exponent*> {
    MATHICGB_ASSERT(mPimpl->ideal.get() != 0);
    MATHICGB_ASSERT(poly < mPimpl->ideal->size());
    const auto& p = *mPimpl->ideal->getPoly(poly);

    MATHICGB_ASSERT(term < p.nTerms());
    return std::make_pair(
      p.coefficientAt(term),
      p.monomialAt(term).unsafeGetRepresentation() + 1
    );
  }
}

// ** Implementation of function mgbi::internalComputeGroebnerBasis
namespace mgbi {
  void internalComputeGroebnerBasis(
    GroebnerInputIdealStream& inputWhichWillBeCleared,
    IdealAdapter& output
  ) {
    /// @todo: make a scheme where the output Groebner basis is freed
    /// polynomial-by-polynomial as data is transferred to out. Also
    /// make it so that ideal is not copied.

    auto&& ideal = PimplOf()(inputWhichWillBeCleared).ideal;
    auto&& conf = inputWhichWillBeCleared.configuration();
    auto&& ring = ideal.ring();
    const auto varCount = ring.getNumVars();
    MATHICGB_ASSERT(PimplOf()(conf).debugAssertValid());

    // Make reducer
    typedef GroebnerConfiguration GConf;
    Reducer::ReducerType reducerType;
    switch (conf.reducer()) {
    case GConf::ClassicReducer:
      reducerType = Reducer::Reducer_BjarkeGeo;
      break;

    default:
    case GConf::DefaultReducer:
    case GConf::MatrixReducer:
      reducerType = Reducer::Reducer_F4_New;
      break;
    }
    const auto reducer = Reducer::makeReducer(reducerType, ring);

    // Set up and configure algorithm
    BuchbergerAlg alg(ideal, 4, *reducer, 2, true, 0);
    alg.setReducerMemoryQuantum(100 * 1024);
    alg.setUseAutoTopReduction(true);
    alg.setUseAutoTailReduction(false);

    // Compute Groebner basis
    alg.computeGrobnerBasis();
    PimplOf()(output).ideal = alg.basis().toIdealAndRetireAll();
  }
}
