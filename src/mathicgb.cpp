#include "mathicgb/stdinc.h"
#include "mathicgb.h"

#include "mathicgb/Basis.hpp"
#include "mathicgb/PolyRing.hpp"
#include "mathicgb/Poly.hpp"
#include "mathicgb/Reducer.hpp"
#include "mathicgb/BuchbergerAlg.hpp"
#include "mathicgb/mtbb.hpp"
#include "mathicgb/LogDomainSet.hpp"
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

namespace mgb {
  bool logTime(const char* logName, double& time) {
    auto log = LogDomainSet::singleton().logDomain(logName);
    if (log == 0 || !log->enabled())
      return false;
    time = log->loggedSecondsReal();
    return true;
  }

  bool logNumber(const char* logName, double& number) {
    auto log = LogDomainSet::singleton().logDomain(logName);
    if (log == 0 || !log->enabled())
      return false;
    number = static_cast<double>(log->count());
    return true;
  }
}

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
      mBaseOrder(ReverseLexicographicBaseOrder),
      mGradings(varCount, 1),
      mReducer(DefaultReducer),
      mMaxSPairGroupSize(0),
      mMaxThreadCount(0),
      mLogging(),
      mCallbackData(0),
      mCallback(0)
#ifdef MATHICGB_DEBUG
      , mHasBeenDestroyed(false)
#endif
    {
    }

    ~Pimpl() {
      MATHICGB_ASSERT(debugAssertValid());
      MATHICGB_IF_DEBUG(mHasBeenDestroyed = true;)
    }

    static bool baseOrderValid(const BaseOrder order) {
      return
        order == LexicographicBaseOrder ||
        order == ReverseLexicographicBaseOrder;
    }

    static bool reducerValid(const Reducer reducer) {
      return
        reducer == DefaultReducer ||
        reducer == ClassicReducer ||
        reducer == MatrixReducer;
    }

    bool debugAssertValid() const {
#ifdef MATHICGB_DEBUG
      MATHICGB_ASSERT(this != 0);
      MATHICGB_ASSERT(baseOrderValid(mBaseOrder));
      MATHICGB_ASSERT(reducerValid(mReducer));
      MATHICGB_ASSERT(mModulus != 0);
      MATHICGB_ASSERT(mCallback != 0 || mCallbackData == 0);
      MATHICGB_ASSERT_NO_ASSUME(!mHasBeenDestroyed);
#endif
      return true;
    }

    const Coefficient mModulus;
    const VarIndex mVarCount;
    BaseOrder mBaseOrder;
    std::vector<Exponent> mGradings;
    Reducer mReducer;
    size_t mMaxSPairGroupSize;
    size_t mMaxThreadCount;
    std::string mLogging;
    void* mCallbackData;
    Callback::Action (*mCallback) (void*);
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

  void GroebnerConfiguration::setMonomialOrderInternal(
    MonomialOrderData order
  ) {
    MATHICGB_ASSERT(Pimpl::baseOrderValid(order.baseOrder));
    MATHICGB_ASSERT(varCount() == 0 || order.gradingsSize % varCount() == 0);
    
    // Currently only reverse lex supported. TODO: make it work
    MATHICGB_ASSERT(order.baseOrder == ReverseLexicographicBaseOrder);

    mPimpl->mBaseOrder = order.baseOrder;
    mPimpl->mGradings.assign
      (order.gradings, order.gradings + order.gradingsSize);
  }

  auto GroebnerConfiguration::monomialOrderInternal() const ->
    MonomialOrderData
  {
    const MonomialOrderData data = {
      mPimpl->mBaseOrder,
      mPimpl->mGradings.data(),
      mPimpl->mGradings.size()
    };
    return data;
  }

  void GroebnerConfiguration::setCallbackInternal(
    void* data,
    Callback::Action (*func) (void*)
  ) {
    MATHICGB_ASSERT(func != 0 || data == 0);
    mPimpl->mCallbackData = data;
    mPimpl->mCallback = func;
  }

  void* GroebnerConfiguration::callbackDataInternal() const {
    return mPimpl->mCallbackData;
  }

  void GroebnerConfiguration::setReducer(Reducer reducer) {
    MATHICGB_ASSERT(Pimpl::reducerValid(reducer));
    mPimpl->mReducer = reducer;
  }

  auto GroebnerConfiguration::reducer() const -> Reducer {
    return mPimpl->mReducer;
  }

  void GroebnerConfiguration::setMaxSPairGroupSize(size_t size) {
    mPimpl->mMaxSPairGroupSize = size;
  }

  size_t GroebnerConfiguration::maxSPairGroupSize() const {
    return mPimpl->mMaxSPairGroupSize;
  }

  void GroebnerConfiguration::setMaxThreadCount(size_t maxThreadCount) {
    mPimpl->mMaxThreadCount = maxThreadCount;
  }

  size_t GroebnerConfiguration::maxThreadCount() const {
    return mPimpl->mMaxThreadCount;
  }

  void GroebnerConfiguration::setLogging(const char* logging) {
    if (logging == 0)
      mPimpl->mLogging.clear();
    else
      mPimpl->mLogging = logging;
  }

  const char* GroebnerConfiguration::logging() const {
    return mPimpl->mLogging.c_str();
  }
}

// ** Implementation of class GroebnerInputIdealStream
namespace mgb {
  struct GroebnerInputIdealStream::Pimpl {
    Pimpl(const GroebnerConfiguration& conf):
      // @todo: varCount should not be int. Fix PolyRing constructor,
      // then remove this static_cast.
      ring(
        conf.modulus(),
        static_cast<int>(conf.varCount()),
        conf.monomialOrder().first ==
          GroebnerConfiguration::LexicographicBaseOrder,
        conf.monomialOrder().second
      ),
      basis(ring),
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
    Basis basis;
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
    MATHICGB_ASSERT(mPimpl->basis.empty());

    MATHICGB_ASSERT(debugAssertValid());
  }

  void GroebnerInputIdealStream::idealBegin(size_t polyCount) {
    MATHICGB_ASSERT(debugAssertValid());
    MATHICGB_IF_DEBUG(mPimpl->checker.idealBegin(polyCount));
    MATHICGB_ASSERT(mPimpl->poly.isZero());
    MATHICGB_ASSERT(mPimpl->basis.empty());

    mPimpl->basis.reserve(polyCount);

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
    mPimpl->basis.insert(std::move(poly));
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
    MATHICGB_ASSERT(&mPimpl->basis.ring() == &mPimpl->ring);
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
    std::unique_ptr<Basis> basis;
  };

  IdealAdapter::IdealAdapter():
    mPimpl(new Pimpl()) {
  }

  IdealAdapter::~IdealAdapter() {
    MATHICGB_ASSERT(mPimpl != 0);
    delete mPimpl;
  }

  auto IdealAdapter::varCount() const -> VarIndex {
    MATHICGB_ASSERT(mPimpl->basis.get() != 0);
    return mPimpl->basis->ring().getNumVars();
  }

  size_t IdealAdapter::polyCount() const {
    MATHICGB_ASSERT(mPimpl->basis.get() != 0);
    return mPimpl->basis->size();
  }

  size_t IdealAdapter::termCount(PolyIndex poly) const {
    MATHICGB_ASSERT(mPimpl->basis.get() != 0);
    MATHICGB_ASSERT(poly < mPimpl->basis->size());
    return mPimpl->basis->getPoly(poly)->nTerms();
  }

  auto IdealAdapter::term(
    PolyIndex poly,
    TermIndex term
  ) const -> std::pair<Coefficient, const Exponent*> {
    MATHICGB_ASSERT(mPimpl->basis.get() != 0);
    MATHICGB_ASSERT(poly < mPimpl->basis->size());
    const auto& p = *mPimpl->basis->getPoly(poly);

    MATHICGB_ASSERT(term < p.nTerms());
    return std::make_pair(
      p.coefficientAt(term),
      p.monomialAt(term).unsafeGetRepresentation() + 1
    );
  }
}

namespace {
  class CallbackAdapter : public BuchbergerAlg::Callback {
  public:
    typedef mgb::GroebnerConfiguration::Callback::Action Action;

    CallbackAdapter(void* data, Action (*callback) (void*)):
      mData(data),
      mCallback(callback),
      mLastAction(Action::ContinueAction)
    {
      MATHICGB_ASSERT(mCallback != 0 || mData == 0);
    }

    const bool isNull() const {return mCallback == 0;}
    Action lastAction() const {return mLastAction;}

    virtual bool call() {
      if (isNull())
        return true;
      mLastAction = mCallback(mData);
      return mLastAction == Action::ContinueAction;
    }

  private:
    void* const mData;
    Action (* const mCallback) (void*);
    Action mLastAction;
  };
}

// ** Implementation of function mgbi::internalComputeGroebnerBasis
namespace mgbi {
  bool internalComputeGroebnerBasis(
    GroebnerInputIdealStream& inputWhichWillBeCleared,
    IdealAdapter& output
  ) {
    /// @todo: make a scheme where the output Groebner basis is freed
    /// polynomial-by-polynomial as data is transferred to out. Also
    /// make it so that ideal is not copied.

    auto&& basis = PimplOf()(inputWhichWillBeCleared).basis;
    auto&& conf = inputWhichWillBeCleared.configuration();
    auto&& ring = basis.ring();
    const auto varCount = ring.getNumVars();
    MATHICGB_ASSERT(PimplOf()(conf).debugAssertValid());

    // Tell tbb how many threads to use
    const auto maxThreadCount = conf.maxThreadCount();
    const auto tbbMaxThreadCount = maxThreadCount == 0 ?
      mgb::tbb::task_scheduler_init::automatic : maxThreadCount;
    mgb::tbb::task_scheduler_init scheduler(tbbMaxThreadCount);

    // Set up logging
    LogDomainSet::singleton().reset();
    LogDomainSet::singleton().performLogCommands(conf.logging());

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
    CallbackAdapter callback(
      PimplOf()(conf).mCallbackData,
      PimplOf()(conf).mCallback
    );

    // Set up and configure algorithm
    BuchbergerAlg alg(basis, 4, *reducer, 2, true, 0);
    alg.setReducerMemoryQuantum(100 * 1024);
    alg.setUseAutoTopReduction(true);
    alg.setUseAutoTailReduction(false);
    alg.setSPairGroupSize(conf.maxSPairGroupSize());
    if (!callback.isNull())
      alg.setCallback(&callback);

    // Compute Groebner basis
    alg.computeGrobnerBasis();
    typedef mgb::GroebnerConfiguration::Callback::Action Action;
    if (callback.lastAction() != Action::StopWithNoOutputAction) {
      PimplOf()(output).basis = alg.basis().toBasisAndRetireAll();
      return true;
    } else
      return false;
  }
}
