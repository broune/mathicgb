#ifndef MATHICGB_MATHICGB_GUARD
#define MATHICGB_MATHICGB_GUARD

#include <ostream>

// The main function in this file is computeGroebnerBasis. See the comment
// preceding that function for an example of how to use this library
// interface.

/// Code in this namespace is not part of the public interface of MathicGB.
/// If you take a dependency on code in this namespace, except indirectly
/// through using code in the mgb namespace, then you are not using the
/// library interface as intended.
namespace mgbi { // Not part of the public interface of MathicGB
  class PimplOf;
}

/// The classes and functions in this namespace make up the public interface
/// of MathicGB. You should not have to update your code when upgrading to a
/// newer minor revision of MathicGB if you only use the public interface.
namespace mgb { // Part of the public interface of MathicGB
  /// Use this class to describe a configuration of a Groebner basis algorithm
  /// that you want to run.
  ///
  /// @todo: expose more of the available functionality.
  class GroebnerConfiguration {
  public:
    typedef unsigned int Coefficient;
    typedef size_t VarIndex;
    typedef int Exponent;

    GroebnerConfiguration(Coefficient modulus, VarIndex varCount);
    GroebnerConfiguration(const GroebnerConfiguration& conf);
    ~GroebnerConfiguration();

    Coefficient modulus() const;
    VarIndex varCount() const;

    enum Reducer {
      DefaultReducer = 0,
      ClassicReducer = 1, /// the classic polynomial division algorithm
      MatrixReducer = 2 /// use linear algebra as in F4
    };
    void setReducer(Reducer reducer);
    Reducer reducer() const;

  private:
    friend class mgbi::PimplOf;

    struct Pimpl;
    Pimpl* mPimpl;
  };

  /// After making a configuration, use this class to communicate the input
  /// basis you want want to run a Groebner basis algorithm on.
  class GroebnerInputIdealStream {
  public:
    GroebnerInputIdealStream(const GroebnerConfiguration& conf);
    ~GroebnerInputIdealStream();

    typedef GroebnerConfiguration::Coefficient Coefficient;
    typedef GroebnerConfiguration::VarIndex VarIndex;
    typedef GroebnerConfiguration::Exponent Exponent;

    const GroebnerConfiguration& configuration();
    Coefficient modulus() const;
    VarIndex varCount() const;

    void idealBegin(size_t polyCount);
    void appendPolynomialBegin(size_t termCount);
    void appendTermBegin();

    /// The sequence of indices appended to a term must be in strictly
    /// ascending order.
    void appendExponent(VarIndex index, Exponent exponent);
    void appendTermDone(Coefficient coefficient);
    void appendPolynomialDone();
    void idealDone();

  private:
    bool debugAssertValid() const;
    Exponent* const mExponents;

    struct Pimpl;
    friend class mgbi::PimplOf;
    Pimpl* const mPimpl;
  };

  /// After making a configuration and an ideal, use this function to compute
  /// a Groebner basis. The output basis is constructed on output, which must
  /// resemble GroebnerInputIdealStream by having the following functions.
  ///
  ///   - modulus() const;
  ///   - varCount() const;
  ///   - idealBegin(size_t polyCount);
  ///   - void appendPolynomialBegin(size_t termCount);
  ///   - void appendTermBegin();
  ///   - void appendExponent(VarIndex index, Exponent exponent);
  ///   - void appendTermDone(Coefficient coefficient);
  ///   - void appendPolynomialDone();
  ///   - void idealDone();
  ///
  /// ** Example
  ///
  /// This example uses a default configuration, constructs an ideal and then
  /// outputs the Groebner basis to a NullIdealStream which does not do
  /// anything with the output. However, we wrap the NullIdealStream in a
  /// IdealStreamLog, which prints out the method calls done on the stream
  /// to std::cerr. We also wrap the GroebnerInputIdealStream in a
  /// IdealStreamChecker which checks that we are correctly following the
  /// protocol of GroebnerInputIdealStream - that's only recommended when
  /// debugging as it is slow.
  ///
  /// GroebnerConfiguration configuration(101, 4); // mod 101, 4 variables
  /// GroebnerInputIdealStream input(configuration);
  /// IdealStreamChecker<GroebnerInputIdealStream> checked(input);
  /// checked.idealBegin(2); // describe ideal with 2 basis elements
  ///   checked.appendPolynomial(2); // describe generator with 2 terms
  ///     checked.appendTermBegin();
  ///       checked.appendExponent(0, 40); // x0^40
  ///     checked.appendTermDone(3); // 3 * x0^40
  ///     checked.appendTermBegin();
  ///       checked.appendExponent(1, 5); // x1^5
  ///       checked.appendExponent(2, 7); // x2^7
  ///     checked.appendTermDone(11); // 11 * x1^5 * x2^7
  ///   checked.appendPolynomialDone(); // 3 * x0^40 + 11 * x1^5 * x2^7
  ///   checked.appendPolynomialBegin(1);
  ///     checked.appendTermBegin();
  ///     checked.appendTermDone(13); // 13
  ///   checked.appendPolynomialDone(); // 13
  /// checked.idealDone(); // the generators are 3*x0^40 + 11*x1^5*x2^7 and 13
  /// NullIdealStream nullStream;
  /// IdealStreamLog<NullIdealStream> logStream(nullStream, std::cerr);
  /// computeGroebnerBasis(input, logStream);
  ///
  /// The ideal constructed on the passed-in GroebnerInputIdealStream will
  /// be cleared. If you need it again for something else, you will have
  /// to re-construct it.
  template<class OutputStream>
  void computeGroebnerBasis(
    GroebnerInputIdealStream& inputWhichWillBeCleared,
    OutputStream& output
  );

  class NullIdealStream;

  /// Passes on all method calls to an inner ideal stream while printing out
  /// what methods get called to an std::ostream.
  template<class Stream = NullIdealStream>
  class IdealStreamLog {
  public:
    typedef GroebnerConfiguration::Coefficient Coefficient;
    typedef GroebnerConfiguration::VarIndex VarIndex;
    typedef GroebnerConfiguration::Exponent Exponent;

    /// All calls are written to log and then passed on to stream.
    IdealStreamLog(std::ostream& log, Stream& stream);

    /// All calls are written to log.
    IdealStreamLog(std::ostream& log, Coefficient modulus, VarIndex varCount);

    ~IdealStreamLog();

    Coefficient modulus() const;
    VarIndex varCount() const;

    void idealBegin(size_t polyCount);
    void appendPolynomialBegin(size_t termCount);
    void appendTermBegin();
    void appendExponent(VarIndex index, Exponent exponent);
    void appendTermDone(Coefficient coefficient);
    void appendPolynomialDone();
    void idealDone();

  private:
    const Coefficient mModulus;
    const VarIndex mVarCount;
    Stream* const mStream;
    std::ostream& mLog;
  };

  /// An ideal stream that simply ignores all the method calls it receives.
  /// Can be handy in combination with IdealStreamLog or as a temporary
  /// stand-in if you have not yet written your own ideal stream.
  class NullIdealStream {
  public:
    typedef GroebnerConfiguration::Coefficient Coefficient;
    typedef GroebnerConfiguration::VarIndex VarIndex;
    typedef GroebnerConfiguration::Exponent Exponent;

    NullIdealStream(Coefficient modulus, VarIndex varCount);

    Coefficient modulus() const {return mModulus;}
    VarIndex varCount() const {return mVarCount;}

    void idealBegin(size_t polyCount) {}
    void appendPolynomialBegin(size_t termCount) {}
    void appendTermBegin() {}
    void appendExponent(VarIndex index, Exponent exponent) {}
    void appendTermDone(Coefficient coefficient) {}
    void appendPolynomialDone() {}
    void idealDone() {}

  private:
    const Coefficient mModulus;
    const VarIndex mVarCount;
  };

  template<class OutputStream>
  void streamSimpleIdeal(OutputStream& output);
}

namespace mgbi { // Not part of the public interface of MathicGB
  using namespace mgb;

  /// Class that checks to see if the protocol for an ideal stream is followed.
  /// A method will return false if that is not the case. This class is not
  /// itself an ideal stream - it is intended to be used with ideal streams for
  /// debugging.
  class StreamStateChecker {
  public:
    typedef GroebnerConfiguration::Coefficient Coefficient;
    typedef GroebnerConfiguration::VarIndex VarIndex;
    typedef GroebnerConfiguration::Exponent Exponent;

    StreamStateChecker(const Coefficient modulus, const VarIndex varCount);
    ~StreamStateChecker();

    void idealBegin(size_t polyCount);
    void appendPolynomialBegin(size_t termCount);
    void appendTermBegin();
    void appendExponent(VarIndex index, Exponent exponent);
    void appendTermDone(Coefficient coefficient);
    void appendPolynomialDone();
    void idealDone();

    bool hasIdeal() const;

  private:
    struct Pimpl;
    Pimpl* const mPimpl;
  };
}

namespace mgb { // Part of the public interface of MathicGB

  /// Use this class to check that you are following the correct protocol
  /// for calling methods on an ideal stream. This has significant overhead
  /// so it is not recommended for production use. If you have built the
  /// MathicGB library in debug mode then this is already automatically
  /// used for GroebnerInputIdealStream.
  template<class Stream>
  class IdealStreamChecker {
  public:
    typedef GroebnerConfiguration::Coefficient Coefficient;
    typedef GroebnerConfiguration::VarIndex VarIndex;
    typedef GroebnerConfiguration::Exponent Exponent;

    IdealStreamChecker(Stream& stream);

    Coefficient modulus() const;
    VarIndex varCount() const;

    void idealBegin(size_t polyCount);
    void appendPolynomialBegin(size_t termCount);
    void appendTermBegin();
    void appendExponent(VarIndex index, Exponent exponent);
    void appendTermDone(Coefficient coefficient);
    void appendPolynomialDone();
    void idealDone();

  private:
    Stream& mStream;
    ::mgbi::StreamStateChecker mChecker;
  };
}

// ********************************************************
// Nothing below this line is part of the public interface of MathicGB.

#ifdef MATHICGB_DEBUG
#include <cassert>
#endif

namespace mgb {
  // ** Functions

  // This method is made inline to avoid the overhead from calling a function
  // for every exponent. This is also why mExponents is not inside the pimpl -
  // otherwise we couldn't access it fro here. That then explains why mExponents
  // is a raw pointer instead of a std::vector - the compiler for the caller
  // and the library must agree on the memory layout of the object and that is
  // less likely to introduce problems for a raw pointer than for a
  // std::vector. In particular, doing it this way allows the library and
  // the caller to use different implementations of the STL.
  inline void GroebnerInputIdealStream::appendExponent(
    const VarIndex index,
    Exponent exponent
  ) {
#ifdef MATHICGB_DEBUG
    assert(index < varCount());
#endif
    mExponents[index] = exponent;
  }


  // ** Implementation of the class IdealStreamLog
  // This class has to be inline as it is a template.

  template<class Stream> 
  IdealStreamLog<Stream>::IdealStreamLog(std::ostream& log, Stream& stream):
    mModulus(stream.modulus()),
    mVarCount(stream.varCount()),
    mStream(&stream),
    mLog(log)
  {
    mLog << "IdealStreamLog s(stream, log); // modulus=" << mModulus
      << ", varCount=" << mVarCount << '\n';
  }

  template<class Stream> 
  IdealStreamLog<Stream>::IdealStreamLog(
    std::ostream& log,
    Coefficient modulus,
    VarIndex varCount
  ):
    mModulus(modulus),
    mVarCount(varCount),
    mStream(0),
    mLog(log)
  {
    mLog << "IdealStreamLog s(stream, " << mModulus << ", " << mVarCount << ");\n";
  }

  template<class Stream> 
  IdealStreamLog<Stream>::~IdealStreamLog() {
    mLog << "// s.~IdealStreamLog();\n";
  }

  template<class Stream> 
  typename IdealStreamLog<Stream>::Coefficient
  IdealStreamLog<Stream>::modulus() const {
    return mModulus;
  }

  template<class Stream> 
  typename IdealStreamLog<Stream>::VarIndex
  IdealStreamLog<Stream>::varCount() const {
    return mVarCount;
  }

  template<class Stream> 
  void IdealStreamLog<Stream>::idealBegin(size_t polyCount) {
    mLog << "s.idealBegin(" << polyCount << "); // polyCount\n";
    if (mStream != 0)
      mStream->idealBegin(polyCount);
  }

  template<class Stream> 
  void IdealStreamLog<Stream>::appendPolynomialBegin(size_t termCount) {
    mLog << "s.appendPolynomialBegin(" << termCount << ");\n";
    if (mStream != 0)
      mStream->appendPolynomialBegin(termCount);
  }

  template<class Stream> 
  void IdealStreamLog<Stream>::appendTermBegin() {
    mLog << "s.appendTermBegin();\n";
    if (mStream != 0)
      mStream->appendTermBegin();
  }

  template<class Stream> 
  void IdealStreamLog<Stream>::appendExponent(VarIndex index, Exponent exponent) {
    mLog << "s.appendExponent(" << index << ", " << exponent <<
      "); // index, exponent\n";
    if (mStream != 0)
      mStream->appendExponent(index, exponent);
  }

  template<class Stream> 
  void IdealStreamLog<Stream>::appendTermDone(Coefficient coefficient) {
    mLog << "s.appendTermDone(" << coefficient << "); // coefficient\n";
    if (mStream != 0)
      mStream->appendTermDone(coefficient);
  }

  template<class Stream> 
  void IdealStreamLog<Stream>::appendPolynomialDone() {
    mLog << "s.appendPolynomialDone();\n";
    if (mStream != 0)
      mStream->appendPolynomialDone();
  }

  template<class Stream> 
  void IdealStreamLog<Stream>::idealDone() {
    mLog << "s.idealDone();\n";
    if (mStream != 0)
      mStream->idealDone();
  }


  // ** Implementation of the class IdealStreamChecker
  // This class has to be inline as it is a template.

  template<class Stream>
  IdealStreamChecker<Stream>::IdealStreamChecker(Stream& stream):
    mStream(stream),
    mChecker(stream.modulus(), stream.varCount())
  {}

  template<class Stream>
  typename IdealStreamChecker<Stream>::Coefficient
  IdealStreamChecker<Stream>::modulus() const {
    return mStream.modulus();
  }

  template<class Stream>
  typename IdealStreamChecker<Stream>::VarIndex
  IdealStreamChecker<Stream>::varCount() const {
    return mStream.varCount();
  }

  template<class Stream>
  void IdealStreamChecker<Stream>::idealBegin(size_t polyCount) {
    mChecker.idealBegin(polyCount);
    mStream.idealBegin(polyCount);
  }

  template<class Stream>
  void IdealStreamChecker<Stream>::appendPolynomialBegin(size_t termCount) {
    mChecker.appendPolynomialBegin(termCount);
    mStream.appendPolynomialBegin(termCount);
  }

  template<class Stream>
  void IdealStreamChecker<Stream>::appendTermBegin() {
    mChecker.appendTermBegin();
    mStream.appendTermBegin();
  }

  template<class Stream>
  void IdealStreamChecker<Stream>::appendExponent(VarIndex index, Exponent exponent) {
    mChecker.appendExponent(index, exponent);
    mStream.appendExponent(index, exponent);
  }

  template<class Stream>
  void IdealStreamChecker<Stream>::appendTermDone(Coefficient coefficient) {
    mChecker.appendTermDone(coefficient);
    mStream.appendTermDone(coefficient);
  }

  template<class Stream>
  void IdealStreamChecker<Stream>::appendPolynomialDone() {
    mChecker.appendPolynomialDone();
    mStream.appendPolynomialDone();
  }

  template<class Stream>
  void IdealStreamChecker<Stream>::idealDone() {
    mChecker.idealDone();
    mStream.idealDone();
  }

  // ** Implementation of the class NullIdealStream
  // This class has to be inline as it is a template.
  inline NullIdealStream::NullIdealStream(
    Coefficient modulus,
    VarIndex varCount
  ):
    mModulus(modulus), mVarCount(varCount) {}
}

namespace mgbi {
  /// Used to read an internal MathicGB ideal without exposing the type of
  /// the ideal.
  class IdealAdapter {
  public:
    typedef GroebnerConfiguration::Coefficient Coefficient;
    typedef GroebnerConfiguration::VarIndex VarIndex;
    typedef GroebnerConfiguration::Exponent Exponent;
    typedef size_t PolyIndex;
    typedef size_t TermIndex;

    IdealAdapter();
    ~IdealAdapter();

    VarIndex varCount() const;
    size_t polyCount() const;
    size_t termCount(PolyIndex poly) const;
    std::pair<Coefficient, const Exponent*> term
      (PolyIndex poly, TermIndex term) const;

  private:
    friend class ::mgbi::PimplOf;
    struct Pimpl;
    Pimpl* mPimpl;
  };

  void internalComputeGroebnerBasis(
    GroebnerInputIdealStream& inputWhichWillBeCleared,
    IdealAdapter& output
  );
}

namespace mgb {
  template<class OutputStream>
  void computeGroebnerBasis(
    GroebnerInputIdealStream& inputWhichWillBeCleared,
    OutputStream& output
  ) {
    mgbi::IdealAdapter ideal;
    mgbi::internalComputeGroebnerBasis(inputWhichWillBeCleared, ideal);

    const auto varCount = ideal.varCount();
    const auto polyCount = ideal.polyCount();
    output.idealBegin(polyCount);
    for (size_t polyIndex = 0; polyIndex < polyCount; ++polyIndex) {
      const auto termCount = ideal.termCount(polyIndex);
      output.appendPolynomialBegin(termCount);
      for (size_t termIndex = 0; termIndex < termCount; ++termIndex) {
        output.appendTermBegin();
        const auto term = ideal.term(polyIndex, termIndex);
        for (size_t var = 0; var < varCount; ++var)
          output.appendExponent(var, term.second[var]);
        output.appendTermDone(term.first);
      }
      output.appendPolynomialDone();
    }
    output.idealDone();
  }
}

#endif
