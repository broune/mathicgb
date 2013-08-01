#include "mathicgb/stdinc.h"

#include <gtest/gtest.h>
#include "mathicgb.h"

namespace {
  template<class Stream>
  void makeBasis(Stream& s) {
    s.idealBegin(2);
      s.appendPolynomialBegin(2); // x^2 - y
        s.appendTermBegin();
          s.appendExponent(0,2);
        s.appendTermDone(1);
        s.appendTermBegin();
          s.appendExponent(1,1);
        s.appendTermDone(s.modulus() - 1);
      s.appendPolynomialDone();
      s.appendPolynomialBegin(2); // x^3-z
        s.appendTermBegin();
          s.appendExponent(0,3);
        s.appendTermDone(1);
        s.appendTermBegin();
          s.appendExponent(2,1);
        s.appendTermDone(s.modulus() - 1);
      s.appendPolynomialDone();
    s.idealDone();
  }

  template<class Stream>
  void makeGroebnerBasis(Stream& s) {
    s.idealBegin(3);
      s.appendPolynomialBegin(2); // x^2 - y
        s.appendTermBegin();
          s.appendExponent(0, 2);
          s.appendExponent(1, 0);
          s.appendExponent(2, 0);
        s.appendTermDone(1);
        s.appendTermBegin();
          s.appendExponent(0, 0);
          s.appendExponent(1, 1);
          s.appendExponent(2, 0);
        s.appendTermDone(s.modulus() - 1);
      s.appendPolynomialDone();
      s.appendPolynomialBegin(2); // xy - z
        s.appendTermBegin();
          s.appendExponent(0, 1);
          s.appendExponent(1, 1);
          s.appendExponent(2, 0);
        s.appendTermDone(1);
        s.appendTermBegin();
          s.appendExponent(0, 0);
          s.appendExponent(1, 0);
          s.appendExponent(2, 1);
        s.appendTermDone(s.modulus() - 1);
      s.appendPolynomialDone();
      s.appendPolynomialBegin(2); // y^2 - xz
        s.appendTermBegin();
          s.appendExponent(0, 0);
          s.appendExponent(1, 2);
          s.appendExponent(2, 0);
        s.appendTermDone(1);
        s.appendTermBegin();
          s.appendExponent(0, 1);
          s.appendExponent(1, 0);
          s.appendExponent(2, 1);
        s.appendTermDone(s.modulus() - 1);
      s.appendPolynomialDone();
    s.idealDone();
  }
}

TEST(MathicGBLib, NullIdealStream) {
  {
    mgb::NullIdealStream stream(2, 3);
    ASSERT_EQ(2, stream.modulus());
    ASSERT_EQ(3, stream.varCount());
    makeBasis(stream);
  }

  {
    mgb::NullIdealStream stream(101, 0);
    ASSERT_EQ(101, stream.modulus());
    ASSERT_EQ(0, stream.varCount());
  }
}

TEST(MathicGBLib, IdealStreamLog) {
  {
    const char* const idealStr = 
      "s.idealBegin(2); // polyCount\n"
      "s.appendPolynomialBegin(2);\n"
      "s.appendTermBegin();\n"
      "s.appendExponent(0, 2); // index, exponent\n"
      "s.appendTermDone(1); // coefficient\n"
      "s.appendTermBegin();\n"
      "s.appendExponent(1, 1); // index, exponent\n"
      "s.appendTermDone(6); // coefficient\n"
      "s.appendPolynomialDone();\n"
      "s.appendPolynomialBegin(2);\n"
      "s.appendTermBegin();\n"
      "s.appendExponent(0, 3); // index, exponent\n"
      "s.appendTermDone(1); // coefficient\n"
      "s.appendTermBegin();\n"
      "s.appendExponent(2, 1); // index, exponent\n"
      "s.appendTermDone(6); // coefficient\n"
      "s.appendPolynomialDone();\n"
      "s.idealDone();\n";

    std::ostringstream out1;
    mgb::IdealStreamLog<> stream1(out1, 7, 3);

    mgb::IdealStreamChecker<decltype(stream1)> checker(stream1);

    std::ostringstream out2;
    mgb::IdealStreamLog<decltype(checker)> stream2(out2, checker);

    std::ostringstream out3;
    mgb::IdealStreamLog<decltype(stream2)> stream3(out3, stream2);

    ASSERT_EQ(7, stream1.modulus());
    ASSERT_EQ(3, stream1.varCount());
    ASSERT_EQ(7, checker.modulus());
    ASSERT_EQ(3, checker.varCount());
    ASSERT_EQ(7, stream2.modulus());
    ASSERT_EQ(3, stream2.varCount());
    ASSERT_EQ(7, stream3.modulus());
    ASSERT_EQ(3, stream3.varCount());

    makeBasis(stream3);
    const auto str1 = std::string(
      "IdealStreamLog s(stream, 7, 3);\n"
    ) + idealStr;
    ASSERT_EQ(str1, out1.str()) << "Displayed expected:\n" << out1.str();

    const auto str2 = std::string(
      "IdealStreamLog s(stream, log); // modulus=7, varCount=3\n"
    ) + idealStr;
    ASSERT_EQ(str2, out2.str()) << "Displayed expected:\n" << out2.str();
    ASSERT_EQ(str2, out3.str()) << "Displayed expected:\n" << out3.str();
  }

  // The ideal <> in no variables
  {
    std::ostringstream out;
    mgb::IdealStreamLog<> stream(out, 101, 0);
    ASSERT_EQ(101, stream.modulus());
    ASSERT_EQ(0, stream.varCount());
    stream.idealBegin(0);
    stream.idealDone();
  }

  // The ideal <1, 0> in no variables
  {
    std::ostringstream out;
    mgb::IdealStreamLog<> stream(out, 101, 0);
    ASSERT_EQ(101, stream.modulus());
    ASSERT_EQ(0, stream.varCount());
    stream.idealBegin(2);
      stream.appendPolynomialBegin(0); // 1
        stream.appendTermBegin();
        stream.appendTermDone(1);
      stream.appendPolynomialDone();
      stream.appendPolynomialBegin(0); // 0
      stream.appendPolynomialDone();
    stream.idealDone();
  }
}

TEST(MathicGBLib, ZeroIdealGB) {
  mgb::GroebnerConfiguration configuration(2, 0);
  mgb::GroebnerInputIdealStream input(configuration);
  std::ostringstream out;
  mgb::IdealStreamLog<> logStream(out, 2, 0);

  input.idealBegin(0);
  input.idealDone();
  mgb::computeGroebnerBasis(input, logStream);

  const auto msg =
    "IdealStreamLog s(stream, 2, 0);\n"
    "s.idealBegin(0); // polyCount\n"
    "s.idealDone();\n";
  EXPECT_EQ(msg, out.str());
}

TEST(MathicGBLib, OneIdealGB) {
  mgb::GroebnerConfiguration configuration(2, 0);
  mgb::GroebnerInputIdealStream input(configuration);
  std::ostringstream out;
  mgb::IdealStreamLog<> logStream(out, 2, 0);

  input.idealBegin(1);
  input.appendPolynomialBegin(1);
  input.appendTermBegin();
  input.appendTermDone(1);
  input.appendPolynomialDone();
  input.idealDone();
  mgb::computeGroebnerBasis(input, logStream);

  const auto msg =
     "IdealStreamLog s(stream, 2, 0);\n"
     "s.idealBegin(1); // polyCount\n"
     "s.appendPolynomialBegin(1);\n"
     "s.appendTermBegin();\n"
     "s.appendTermDone(1); // coefficient\n"
     "s.appendPolynomialDone();\n"
     "s.idealDone();\n";
  EXPECT_EQ(msg, out.str());
}

TEST(MathicGBLib, EasyGB) {
  mgb::GroebnerConfiguration configuration(101, 3);
  mgb::GroebnerInputIdealStream input(configuration);
  std::ostringstream computedStr;
  mgb::IdealStreamLog<> computed(computedStr, 101, 3);
  mgb::IdealStreamChecker<decltype(computed)> checked(computed);

  makeBasis(input);
  mgb::computeGroebnerBasis(input, checked);

  std::ostringstream correctStr;
  mgb::IdealStreamLog<> correct(correctStr, 101, 3);
  mgb::IdealStreamChecker<decltype(correct)> correctChecked(correct);
  makeGroebnerBasis(correctChecked);

  EXPECT_EQ(correctStr.str(), computedStr.str())
    << "\nDisplayed expected:\n" << correctStr.str()
    << "\nDisplayed computed:\n" << computedStr.str();
}

TEST(MathicGBLib, EasyReGB) {
  mgb::GroebnerConfiguration configuration(101, 3);
  mgb::GroebnerInputIdealStream input(configuration);
  std::ostringstream computedStr;
  mgb::IdealStreamLog<> computed(computedStr, 101, 3);
  mgb::IdealStreamChecker<decltype(computed)> checked(computed);

  makeGroebnerBasis(input);
  mgb::computeGroebnerBasis(input, checked);

  std::ostringstream correctStr;
  mgb::IdealStreamLog<> correct(correctStr, 101, 3);
  mgb::IdealStreamChecker<decltype(correct)> correctChecked(correct);
  makeGroebnerBasis(correctChecked);

  EXPECT_EQ(correctStr.str(), computedStr.str())
    << "\nDisplayed expected:\n" << correctStr.str()
    << "\nDisplayed computed:\n" << computedStr.str();
}
