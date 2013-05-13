// Copyright 2011 Michael E. Stillman

#include "mathicgb/stdinc.h"
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>

#include "mathicgb/Poly.hpp"
#include "mathicgb/Basis.hpp"
#include "mathicgb/MonTableNaive.hpp"
#include "mathicgb/MonTableKDTree.hpp"
#include "mathicgb/MonTableDivList.hpp"
#include "mathicgb/MTArray.hpp"
#include "mathicgb/io-util.hpp"

#include "mathicgb/MonomialHashTable.hpp"
#include "mathicgb/PolyHashTable.hpp"
#include "mathicgb/PolyHeap.hpp"
#include "mathicgb/PolyGeoBucket.hpp"
#include "mathicgb/SigPolyBasis.hpp"
#include "mathicgb/SignatureGB.hpp"

#include <gtest/gtest.h>

std::string ideal1 =
"32003 4 1 1 1 1 1 \
2 \
ab+3d2+2d \
c3+a-b+1 \
";

std::string ideal2 = " \
32003 4 \
1 1 1 1 1 \
7 \
bc2-ad2 \
abc-d3 \
b3-acd \
a2d2-cd3 \
ac3d-ab2d2 \
a2c2d-b2d3 \
c3d3-b2d4 \
";

TEST(PolyRing, read) {
  std::stringstream o;
  std::string ringinfo = "32003 6\n1 1 1 1 1 1";
  std::unique_ptr<PolyRing> R(ringFromString(ringinfo));
  R->write(o, true);

  EXPECT_EQ("32003 6\nrevlex 1\n 1 1 1 1 1 1\n", o.str());
}

TEST(Poly,readwrite) {
  std::string f1 = "14ce2<72>+13adf<16>";
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  Poly f(*R);
  std::stringstream ifil(f1);
  f.parseDoNotOrder(ifil);
  std::ostringstream o;
  f.display(o,true);
  EXPECT_EQ(o.str(), f1);
}

bool testPolyParse(PolyRing* R, std::string s)
{
  // parse poly, then see if it matches the orig string
  Poly f(*R);
  std::istringstream i(s);
  f.parseDoNotOrder(i);
  std::ostringstream o;
  f.display(o);
  //  std::cout << "orig = " << s << std::endl;
  //  std::cout << "f    = " << o.str() << std::endl;
  return o.str() == s;
}
bool testPolyParse2(PolyRing* R, std::string s, std::string answer)
{
  // parse poly, then see if it matches the orig string
  Poly f(*R);
  std::istringstream i(s);
  f.parseDoNotOrder(i);
  std::ostringstream o;
  f.display(o);
  //  std::cout << "orig = " << s << std::endl;
  //  std::cout << "f    = " << o.str() << std::endl;
  return o.str() == answer;
}

TEST(Poly,parse) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  EXPECT_TRUE(testPolyParse(R.get(), "3a<1>+<0>"));
  EXPECT_TRUE(testPolyParse(R.get(), "3a<1>+13af3<0>+14cde<0>"));
  EXPECT_TRUE(testPolyParse(R.get(), "<1>+13af3<0>+14cde<0>"));
}

bool testMonomialParse(PolyRing* R, std::string s)
{
  Monomial m = stringToMonomial(R, s);
  std::string str2 = monomialToString(R,m);
  return s == str2;
}

TEST(Monomial, parse) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  EXPECT_TRUE(testMonomialParse(R.get(), "ab2d<2>"));
  EXPECT_TRUE(testMonomialParse(R.get(), "ab2d<0>"));
  EXPECT_TRUE(testMonomialParse(R.get(), "<13>"));
  EXPECT_TRUE(testMonomialParse(R.get(), "abcdef<0>"));
  EXPECT_TRUE(testMonomialParse(R.get(), "a10b3d4<0>"));
}

TEST(Monomial,compare) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  
  Monomial mone = stringToMonomial(R.get(), "<0>");
  Monomial mone2 = stringToMonomial(R.get(), "1");
  Monomial m1 = stringToMonomial(R.get(), "ab2<0>");
  Monomial m2 = stringToMonomial(R.get(), "a2b<0>");

  EXPECT_TRUE(R->monomialEQ(mone, mone2));

  //  monomial mone = monomialFromString(R, "0 0");
  //  monomial m1 = monomialFromString(R, "0 2 1 2 0 1");
  //  monomial m2 = monomialFromString(R, "0 2 1 1 0 2");

  bool a = R->monomialLT(m1,m2);
  EXPECT_TRUE(a);
  EXPECT_FALSE(R->monomialEQ(m1,m2));
  EXPECT_EQ(LT, R->monomialCompare(m1,m2));

  a = R->monomialLT(mone,m1);
  EXPECT_TRUE(a);
  EXPECT_FALSE(R->monomialEQ(mone,m1));
  EXPECT_EQ(GT, R->monomialCompare(m1,mone));

  a = R->monomialLT(mone,mone);
  EXPECT_FALSE(a);
  EXPECT_TRUE(R->monomialEQ(mone,mone));
  EXPECT_EQ(EQ, R->monomialCompare(mone,mone));

  Monomial b = stringToMonomial(R.get(), "b<0>");
  Monomial c = stringToMonomial(R.get(), "c<0>");

  a = R->monomialLT(b,c);
  EXPECT_FALSE(a);
}

TEST(Monomial,mult) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial m1 = stringToMonomial(R.get(), "ab2<0>");
  Monomial m2 = stringToMonomial(R.get(), "a2b<0>");
  Monomial m3ans = stringToMonomial(R.get(), "a3b3<0>");

  Monomial m3 = R->allocMonomial();
  R->monomialMult(m1,m2,m3);
  EXPECT_TRUE(R->monomialEQ(m3ans,m3));

  R->freeMonomial(m1);
  R->freeMonomial(m2);
  R->freeMonomial(m3);
  R->freeMonomial(m3ans);
}

TEST(Monomial,multTo) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial m1 = stringToMonomial(R.get(), "ab2<0>");
  Monomial m2 = stringToMonomial(R.get(), "a2b<0>");
  Monomial m3ans = stringToMonomial(R.get(), "a3b3<0>");
  R->monomialMultTo(m1,m2);
  EXPECT_TRUE(R->monomialEQ(m3ans,m1));

  R->freeMonomial(m1);
  R->freeMonomial(m2);
  R->freeMonomial(m3ans);
}

TEST(Monomial, divide) {
  // test of monomialDivide, monomialIsDivisibleBy, monomialQuotientAndMult
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial m1 = stringToMonomial(R.get(), "ab2<0>");
  Monomial m2 = stringToMonomial(R.get(), "a2b<0>");
  Monomial m3ans = stringToMonomial(R.get(), "a3b3<0>");
  Monomial m3 = R->allocMonomial();
  Monomial m1a = R->allocMonomial();

  R->monomialMult(m1,m2,m3);
  EXPECT_TRUE(R->monomialIsDivisibleBy(m3,m2));
  EXPECT_FALSE(R->monomialIsDivisibleBy(m2,m3));
  R->monomialDivide(m3,m2,m1a);
  EXPECT_TRUE(R->monomialEQ(m1,m1a));

  R->freeMonomial(m1);
  R->freeMonomial(m2);
  R->freeMonomial(m3);
  R->freeMonomial(m3ans);
  R->freeMonomial(m1a);
}

TEST(Monomial, monomialQuotientAndMult) {
  // test of monomialQuotientAndMult
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial m1 = stringToMonomial(R.get(), "ab2f2<0>");
  Monomial m2 = stringToMonomial(R.get(), "af<0>");
  Monomial m3 = stringToMonomial(R.get(), "<2>");

  Monomial n = R->allocMonomial();
  Monomial n1 = R->allocMonomial();
  Monomial na = R->allocMonomial();

  R->monomialQuotientAndMult(m1,m2,m3,n);
  R->monomialDivide(m1,m2,n1);  // m1//m2
  R->monomialMult(n1,m3,na); // m1//m2 * m3  should be n

  EXPECT_TRUE(R->monomialEQ(n, na));

  R->freeMonomial(m1);
  R->freeMonomial(m2);
  R->freeMonomial(m3);
  R->freeMonomial(n);
  R->freeMonomial(n1);
  R->freeMonomial(na);
}

void testMonomialOps(const PolyRing* R, std::string s1, std::string s2)
{
  Monomial m1 = stringToMonomial(R, s1);
  Monomial m2 = stringToMonomial(R, s2);
  Monomial m3 = stringToMonomial(R, "abcdef<0>");

 
  Monomial m4 = R->allocMonomial();
  Monomial lcm = R->allocMonomial();
  Monomial m8 = R->allocMonomial();
  Monomial m1a = R->allocMonomial();
  Monomial m2a = R->allocMonomial();
  Monomial m1b = R->allocMonomial();
  Monomial m2b = R->allocMonomial();

  R->monomialMult(m1,m2,m4);
  R->monomialLeastCommonMultiple(m1,m2,lcm);

  // lcm(m1,m2)/m1, lcm(m1,m2)/m2:  relatively prime
  EXPECT_TRUE(R->monomialIsDivisibleBy(lcm, m1));
  EXPECT_TRUE(R->monomialIsDivisibleBy(lcm, m2));
  R->monomialDivide(lcm, m1, m1a);
  R->monomialDivide(lcm, m2, m2a);
  EXPECT_TRUE(R->monomialRelativelyPrime(m1a,m2a));

  EXPECT_TRUE(R->monomialIsDivisibleBy(lcm, m1a));
  EXPECT_TRUE(R->monomialIsDivisibleBy(lcm, m2a));
  R->monomialDivide(lcm, m1a, m1b);
  R->monomialDivide(lcm, m2a, m2b);
  EXPECT_TRUE(R->monomialEQ(m1, m1b));
  EXPECT_TRUE(R->monomialEQ(m2, m2b));
  R->monomialMult(m1a, m2a, m8);

  size_t supp1 = R->monomialSizeOfSupport(m1a);
  size_t supp2 = R->monomialSizeOfSupport(m2a);
  size_t supp = R->monomialSizeOfSupport(m8);
  EXPECT_EQ(supp1+supp2, supp)
    << monomialToString(R,m1) << " " << monomialToString(R,m2) << "\n"
    << monomialToString(R,m1a) << " " << monomialToString(R,m2a) << " "
    << monomialToString(R,m8);
}

TEST(Monomial, ops)
{
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  testMonomialOps(R.get(), "ab2f2<0>", "bc2df3<0>");
  testMonomialOps(R.get(), "ab2f2<0>", "<0>");
  testMonomialOps(R.get(), "<0>", "<0>");
  testMonomialOps(R.get(), "a<0>", "a<0>");
  testMonomialOps(R.get(), "a<0>", "b<0>");
  testMonomialOps(R.get(), "a10b10c10d10e10f10<0>", "b2f5<0>");
}

TEST(Monomial, ei)
{
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial m1 = stringToMonomial(R.get(), "<1>");
  Monomial m1a = R->allocMonomial();
  R->monomialEi(1, m1a);
  EXPECT_TRUE(R->monomialEQ(m1,m1a));


  m1 = stringToMonomial(R.get(), "<0>");
  R->monomialEi(0, m1a);
  EXPECT_TRUE(R->monomialEQ(m1,m1a));

  m1 = stringToMonomial(R.get(), "<1>");
  R->monomialEi(1, m1a);
  EXPECT_TRUE(R->monomialEQ(m1,m1a));

  m1 = stringToMonomial(R.get(), "<10000>");
  R->monomialEi(10000, m1a);
  EXPECT_TRUE(R->monomialEQ(m1,m1a));
}

TEST(Monomial, strict)
{
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial m1 = stringToMonomial(R.get(), "ab2c3d4e");
  Monomial m2 = stringToMonomial(R.get(), "ab2c3d4");
  Monomial m3 = stringToMonomial(R.get(), "ab2c3d4");
  Monomial m4 = stringToMonomial(R.get(), "ab2c3d3e");

  EXPECT_TRUE(R->monomialHasStrictlyLargerExponent(m1,m2,m3));
  EXPECT_FALSE(R->monomialHasStrictlyLargerExponent(m1,m2,m4));
}

TEST(Monomial, divideToNegative)
{
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial m1 = stringToMonomial(R.get(), "ab100<0>");
  Monomial m2 = stringToMonomial(R.get(), "ab2c3d4<0>");
  Monomial m3 = R->allocMonomial();
  Monomial m4 = R->allocMonomial();
  Monomial m5 = R->allocMonomial();
  Monomial mone = stringToMonomial(R.get(), "<0>");

  R->monomialDivideToNegative(m1,m2,m3);
  R->monomialDivideToNegative(m2,m1,m4);
  R->monomialMult(m3,m4,m5);

  EXPECT_TRUE(R->monomialEQ(m5,mone));

  m3 = stringToMonomial(R.get(), "ab2c3d4");
  m4 = stringToMonomial(R.get(), "ab2c3d3e");

  EXPECT_TRUE(R->monomialHasStrictlyLargerExponent(m1,m2,m3));
  EXPECT_FALSE(R->monomialHasStrictlyLargerExponent(m2,m2,m4));
}

TEST(Monomial, findSignature)
{
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));

  Monomial v1 = stringToMonomial(R.get(), "abef");
  Monomial v2 = stringToMonomial(R.get(), "acdf2");
  Monomial u1 = stringToMonomial(R.get(), "f5<13>");
  Monomial t1 = R->allocMonomial();
  Monomial t1ans = stringToMonomial(R.get(), "cdf6<13>");

  R->monomialFindSignature(v1,v2,u1,t1);
  EXPECT_TRUE(R->monomialEQ(t1,t1ans));
}

//#warning "remove this code"
#if 0
bool testMonomialOldParse(PolyRing *R, std::string s)
{
  monomial m = monomialParseFromString(R, s);
  std::string str2 = monomialDisplay(R,m);
  return s == str2;
}

TEST(OldMonomial, parse) {
  PolyRing *R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  EXPECT_TRUE(testMonomialOldParse(R, "ab2d<2>"));
  EXPECT_TRUE(testMonomialOldParse(R, "ab2d<0>"));
  EXPECT_TRUE(testMonomialOldParse(R, "<13>"));
  EXPECT_TRUE(testMonomialOldParse(R, "abcdef<0>"));
  EXPECT_TRUE(testMonomialOldParse(R, "a10b3d4<0>"));
}

TEST(OldMonomial,compare) {
  PolyRing *R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  monomial mone = monomialFromString(R, "0 0");
  monomial m1 = monomialFromString(R, "0 2 1 2 0 1");
  monomial m2 = monomialFromString(R, "0 2 1 1 0 2");

  bool a = R->monomialLT(m1,m2);
  EXPECT_TRUE(a);
  EXPECT_FALSE(R->monomialEQ(m1,m2));
  EXPECT_EQ(LT, R->monomialCompare(m1,m2));

  a = R->monomialLT(mone,m1);
  EXPECT_TRUE(a);
  EXPECT_FALSE(R->monomialEQ(mone,m1));
  EXPECT_EQ(GT, R->monomialCompare(m1,mone));

  a = R->monomialLT(mone,mone);
  EXPECT_FALSE(a);
  EXPECT_TRUE(R->monomialEQ(mone,mone));
  EXPECT_EQ(EQ, R->monomialCompare(mone,mone));

  monomial b = monomialFromString(R, "0 1 1 1");
  monomial c = monomialFromString(R, "0 1 2 1");
  a = R->monomialLT(b,c);
  EXPECT_FALSE(a);
}

TEST(OldMonomial,mult) {
  PolyRing *R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  monomial m1 = monomialFromString(R, "0 2 1 2 0 1"); // ab2
  monomial m2 = monomialFromString(R, "0 2 1 1 0 2"); // a2b
  // std::cout << "m1 is " << monomialToString(R,m1) << std::endl;
  //std::cout << "m2 is " << monomialToString(R,m2) << std::endl;
  monomial m3 = new int[R->maxMonomialSize()];
  R->monomialMult(m1,m2,m3);
  monomial m3ans = monomialFromString(R, "0 2 1 3 0 3"); // a3b3
  //std::cout << "answer is " << monomialToString(R,m3) << std::endl;
  EXPECT_TRUE(R->monomialEQ(m3ans,m3));
}

TEST(OldMonomial,multTo) {
  PolyRing *R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  monomial m1 = monomialFromString(R, "0 3 5 2 1 2 0 1"); // ab2
  monomial m2 = monomialFromString(R, "0 2 1 1 0 2"); // a2b
  R->monomialMultTo(m1,m2);
  monomial m3ans = monomialFromString(R, "0 3 5 2 1 3 0 3"); // a3b3
  EXPECT_TRUE(R->monomialEQ(m3ans,m1));
}

TEST(OldMonomial, divide) {
  // test of monomialDivide, monomialIsDivisibleBy, monomialQuotientAndMult
  PolyRing *R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  monomial m1 = monomialFromString(R, "0 3 5 2 1 2 0 1"); // ab2
  monomial m2 = monomialFromString(R, "0 2 1 1 0 2"); // a2b
  monomial m3 = new int[R->maxMonomialSize()];
  monomial m1a = new int[R->maxMonomialSize()];
  R->monomialMult(m1,m2,m3);
  EXPECT_TRUE(R->monomialIsDivisibleBy(m3,m2));
  EXPECT_FALSE(R->monomialIsDivisibleBy(m2,m3));
  R->monomialDivide(m3,m2,m1a);
  EXPECT_TRUE(R->monomialEQ(m1,m1a));
}

TEST(OldMonomial, monomialQuotientAndMult) {
  // test of monomialQuotientAndMult
  PolyRing *R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  monomial m1 = monomialFromString(R, "0 3 5 2 1 2 0 1"); // ab2f2
  monomial m2 = monomialFromString(R, "0 2 5 1 0 1"); // af
  monomial m3 = monomialFromString(R, "2 0"); // e_2
  monomial n = new int[R->maxMonomialSize()];
  monomial n1 = new int[R->maxMonomialSize()];
  monomial na = new int[R->maxMonomialSize()];
  R->monomialQuotientAndMult(m1,m2,m3,n);
  R->monomialDivide(m1,m2,n1);  // m1//m2
  R->monomialMult(n1,m3,na); // m1//m2 * m3  should be n
  //  std::cout << "n is " << monomialToString(R,n) << std::endl;
  //  std::cout << "na is " << monomialToString(R,na) << std::endl;

  EXPECT_TRUE(R->monomialEQ(n, na));
}
#endif


TEST(Coeff, reciprocal) {
  std::unique_ptr<PolyRing> R(ringFromString("11 6 1\n1 1 1 1 1 1"));
  coefficient vals[10] = {1,2,3,4,5,6,7,8,9,10};
  coefficient ans[10] = {1,6,4,3,9,2,8,7,5,10};
  for (int i=0; i<10; i++)
    {
      coefficient a = vals[i];
      R->coefficientReciprocalTo(a);
      EXPECT_EQ(ans[i], a);
    }
}

TEST(Coeff, reciprocal2) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  for (int i=1; i<32002; i++)
    {
      coefficient a1;
      coefficient a = i;
      R->coefficientReciprocalTo(a);
      R->coefficientDivide(1,i,a1);
      EXPECT_EQ(a,a1);
      R->coefficientMultTo(a,i);
      EXPECT_TRUE(a > 0);
      EXPECT_TRUE(a < 32003);
      EXPECT_TRUE(a == 1);
    }
}

TEST(Coeff, negate) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  for (int i=1; i<32002; i++)
    {
      coefficient a1 = i;
      coefficient a = i;
      R->coefficientNegateTo(a);
      EXPECT_TRUE(a > 0);
      EXPECT_TRUE(a < 32003);
      R->coefficientAddTo(a1,a);
      EXPECT_EQ(0,a1);
    }
}

TEST(Coeff, addone) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  for (int i=0; i<32002; i++)
    {
      coefficient a1 = i;
      coefficient a = i;
      R->coefficientAddOneTo(a);
      EXPECT_TRUE(a > 0);
      EXPECT_TRUE(a < 32003);
      R->coefficientAddTo(a1,1);
      EXPECT_EQ(a,a1);
    }
}

TEST(MTArray,Naive1) {
  // We create a table here
  size_t not_used = 0;
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  std::unique_ptr<MonomialTableArray> M(MonomialTableArray::make(R.get(), 1, 6, false));
  std::string mons[2] = {
    "abc<1>",
    "a2d<1>"
  };
  for (int i=0; i<2; i++)
    {
      monomial m = monomialParseFromString(R.get(), mons[i]);
      M->insert(m,0);
    }
  //  M.display(std::cout);

  // Now we test membership
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "abc4d<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "a2d2<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "abc<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "a2d<1>"),not_used));
  EXPECT_FALSE(M->member(monomialParseFromString(R.get(), "a2d<2>"),not_used));
  EXPECT_FALSE(M->member(monomialParseFromString(R.get(), "ad<1>"),not_used));
}


TEST(MTArray,DivList1) {
  // We create a table here
  size_t not_used = 0;
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  auto M = MonomialTableArray::make(R.get(), 1, 6, false);
  std::string mons[2] = {
    "abc<1>",
    "a2d<1>"
  };
  for (int i=0; i<2; i++)
    {
      monomial m = monomialParseFromString(R.get(), mons[i]);
      M->insert(m,0);
    }
  //  M.display(std::cout);

  // Now we test membership
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "abc4d<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "a2d2<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "abc<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "a2d<1>"),not_used));
  EXPECT_FALSE(M->member(monomialParseFromString(R.get(), "a2d<2>"),not_used));
  EXPECT_FALSE(M->member(monomialParseFromString(R.get(), "ad<1>"),not_used));
}

TEST(MTArray,KDTree1) {
  // We create a table here
  size_t not_used = 0;
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  std::unique_ptr<MonomialTableArray> M(MonomialTableArray::make(R.get(), 2, 6, false));
  std::string mons[2] = {
    "abc<1>",
    "a2d<1>"
  };
  for (int i=0; i<2; i++)
    {
      monomial m = monomialParseFromString(R.get(), mons[i]);
      M->insert(m,0);
    }
  //  M.display(std::cout);

  // Now we test membership
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "abc4d<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "a2d2<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "abc<1>"),not_used));
  EXPECT_TRUE(M->member(monomialParseFromString(R.get(), "a2d<1>"),not_used));
  EXPECT_FALSE(M->member(monomialParseFromString(R.get(), "a2d<2>"),not_used));
  EXPECT_FALSE(M->member(monomialParseFromString(R.get(), "ad<1>"),not_used));
}

//#warning "remove this code"
#if 0
bool test_find_signatures(const PolyRing *R, 
			  const_monomial u1, 
			  const_monomial u2, 
			  const_monomial v1, 
			  const_monomial v2)
{
  monomial g = new int[R->maxMonomialSize()];
  monomial t1 = new int[R->maxMonomialSize()];
  monomial t2 = new int[R->maxMonomialSize()];
  monomial x1 = new int[R->maxMonomialSize()];
  monomial x2 = new int[R->maxMonomialSize()];
  monomial y1 = new int[R->maxMonomialSize()];
  monomial y2 = new int[R->maxMonomialSize()];
  monomial v1v2 = new int[R->maxMonomialSize()];
  monomial x1g = new int[R->maxMonomialSize()];
  monomial p = new int[R->maxMonomialSize()];
  monomial m = new int[R->maxMonomialSize()];
  

  R->monomialFindSignature(v1,v2,u1,t1);
  R->monomialFindSignature(v2,v1,u2,t2);

  R->monomialMult(v1, v2, p);
  R->monomialLeastCommonMultiple(v1, v2, m);
  R->monomialDivide(p, m, g);
  //R->monomialGreatestCommonDivisor(v1, v2, g);

  // check that v1*t1 == v2*t2
  // v1*v2 == g * (v1*t1)
  R->monomialDivide(t1,u1,y1);
  R->monomialDivide(t2,u2,y2); // remove mult by signatures
  R->monomialMult(v1,y1,x1);
  R->monomialMult(v2,y2,x2);
  R->monomialMult(v1,v2,v1v2);
  R->monomialMult(x1,g,x1g);
#if 0
  std::cout << "v1 is " << monomialToString(R,v1) << std::endl;
  std::cout << "v2 is " << monomialToString(R,v2) << std::endl;
  std::cout << "u1 is " << monomialToString(R,u1) << std::endl;
  std::cout << "u2 is " << monomialToString(R,u2) << std::endl;
  std::cout << "t1 is " << monomialToString(R,t1) << std::endl;
  std::cout << "22 is " << monomialToString(R,t2) << std::endl;
  std::cout << "x1 is " << monomialToString(R,x1) << std::endl;
  std::cout << "x2 is " << monomialToString(R,x2) << std::endl;
  std::cout << "y1 is " << monomialToString(R,y1) << std::endl;
  std::cout << "y2 is " << monomialToString(R,y2) << std::endl;
  std::cout << "g is " << monomialToString(R,g) << std::endl;
  std::cout << "v1v2 is " << monomialToString(R,v1v2) << std::endl;
  std::cout << "x1g is " << monomialToString(R,x1g) << std::endl;
#endif
  if (!R->monomialEQ(x1,x2)) return false;
  if (!R->monomialEQ(v1v2,x1g)) return false;

  return true;
}

TEST(OldMonomial, findSignatures) {
  PolyRing *R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  monomial v1 = monomialFromString(R, "0 3  5 2  1 2  0 1"); // ab2f2
  monomial v2 = monomialFromString(R, "0 3  2 3  1 1  0 1"); // abc3
  monomial u1 = monomialFromString(R, "2 3  4 1  1 2  0 1"); // 
  monomial u2 = monomialFromString(R, "3 2  1 1  0 2"); // 
  EXPECT_TRUE(test_find_signatures(R,u1,u2,v1,v2));
}
#endif

TEST(Ideal,readwrite) {
  // This also tests Poly::iterator
  std::unique_ptr<Basis> I = basisParseFromString(ideal1);
  size_t ngens = I->viewGenerators().size();
  EXPECT_TRUE(2 == ngens);

  // now read and write each generator
  for (size_t i=0; i<ngens; i++)
    {
      const Poly *f = I->getPoly(i);
      std::ostringstream o;
      f->display(o,false);
      Poly g(f->ring());
      std::stringstream ifil(o.str());
      g.parse(ifil);
      EXPECT_TRUE(g == *f);
    }
}

TEST(Poly,lead) {
  // This also tests Poly::iterator, Poly::read, Poly::write
  std::unique_ptr<Basis> I = basisParseFromString(ideal1);
  std::unique_ptr<const PolyRing> R(I->getPolyRing());
  monomial lm = stringToMonomial(R.get(), "ab");
  EXPECT_TRUE(R->monomialEQ(lm, I->getPoly(0)->getLeadMonomial()));
  EXPECT_EQ(1, I->getPoly(0)->getLeadCoefficient());
  EXPECT_EQ(0, I->getPoly(0)->getLeadComponent());
  R->freeMonomial(lm);
}

//////////////////////////////
// Test reducer code /////////
//////////////////////////////

std::unique_ptr<Poly> multIdealByPolyReducer(int typ, const Basis& basis, const Poly& g)
{
  const PolyRing& R = basis.ring();
  auto poly = make_unique<Poly>(R);
  std::unique_ptr<Reducer> H = Reducer::makeReducer(static_cast<Reducer::ReducerType>(typ), R);
  for (Poly::const_iterator i = g.begin(); i != g.end(); ++i) {
    monomial mon = R.allocMonomial();
    R.monomialCopy(i.getMonomial(), mon);
    int x = R.monomialGetComponent(mon);
    R.monomialChangeComponent(mon, 0);
    std::unique_ptr<Poly> h(basis.getPoly(x)->copy());
    h->multByTerm(i.getCoefficient(), mon);
    R.monomialSetIdentity(mon);

    size_t ncmps;
    Poly* sum =
      Poly::add(&R, h->begin(), h->end(), poly->begin(), poly->end(), ncmps); 
    poly.reset(sum);
  }
  return poly;
}

void testPolyReducer(
  Reducer::ReducerType reducerType,
  const Basis& basis,
  const std::string& f,
  const std::string& ans
) {
  const PolyRing& ring = *basis.getPolyRing();
  std::unique_ptr<Poly> g = polyParseFromString(&ring, f);
  std::unique_ptr<Poly> h = multIdealByPolyReducer(reducerType, basis, *g);
  if (!h->isZero()) {
    Poly::iterator prev = h->begin();
    Poly::iterator it = prev;
    for (++it; it != h->end(); ++it, ++prev) {
      EXPECT_TRUE(ring.monomialLT(it.getMonomial(), prev.getMonomial()))
        << "Reduced result not in sorted order: " << toString(h.get());
    }
  }

  EXPECT_EQ(ans, toString(h.get())) << "Reducer type " << reducerType;
}

/// @todo: this is no longer a test of a reducer. What to do this this test?
TEST(Reducer, insert) {
  // plan: read an ideal, and another poly
  //  use this last poly to determine what to add to the heap
  // at the end, take the value of the heap, compare to answer

  std::unique_ptr<Basis> I = basisParseFromString(ideal2); // charac is 32003
  for (int typ = 0; typ <= 30; ++typ) {
    Reducer::ReducerType red = Reducer::ReducerType(typ);
    if (static_cast<int>(red) != typ ||
      Reducer::makeReducerNullOnUnknown(red, I->ring()).get() == 0)
      continue;

    testPolyReducer(red, *I, "c2d<0>", "bc4d<0>-ac2d3<0>");
    testPolyReducer(red, *I, "c2d<0>-3abc<1>", "-3a2b2c2<0>+bc4d<0>+3abcd3<0>-ac2d3<0>");
    testPolyReducer(red, *I, "c2d<0>-3abc<1>+a2<3>", "-3a2b2c2<0>+bc4d<0>+a4d2<0>-a2cd3<0>+3abcd3<0>-ac2d3<0>");
    testPolyReducer(red, *I, "c2d<0>-3abc<1>+a2<3>+3ab<4>", "3a2bc3d<0>-3a2b3d2<0>-3a2b2c2<0>+bc4d<0>+a4d2<0>-a2cd3<0>+3abcd3<0>-ac2d3<0>");
    testPolyReducer(red, *I, "bc4d<0>+3bc4d<0>-2bc4d<0>", "2b2c6d<0>-2abc4d3<0>");
    testPolyReducer(red, *I, "a<0>+32002c<1>+<3>", "0");
  }
}

std::string somePolys =
"  bc+bd+af+ef\n\
  ac+cd+bf+ef\n\
  ad+bd+cf+ef\n\
  ab+ac+df+ef\n\
  bd2+cd2+c2f+adf+bdf+cef\n\
  b2d+acd+bcf+d2f+bef+def\n\
  c2d+bd2+b2f+bcf+adf+cdf+bef+def\n\
  a2f+b2f+c2f+d2f+aef+bef+cef+def\n\
  cd2f+d3f+adef+bdef+cdef+d2ef+b2f2+acf2+c2f2+bdf2+cdf2+aef2+cef2+def2\n\
  c3f+b2df+acdf+cd2f+c2ef+bdef\n\
  b3f+b2cf+bc2f+b2df+acdf+bcdf+b2ef+bcef+bdef+cdef\n\
  b2ef2+bcef2+adef2+d2ef2+be2f2+de2f2+abf3+b2f3+bcf3+adf3+cdf3+d2f3\n\
  bde2f2+cde2f2+abdf3+b2df3+bcdf3+ad2f3+cd2f3+d3f3+b2ef3+acef3+bcef3+c2ef3+adef3+d2ef3+ae2f3+be2f3+de2f3+e3f3\n\
  c2e2f3+ade2f3+bde2f3+d2e2f3+ce3f3+de3f3+b2df4+acdf4+bcdf4+c2df4+ad2f4+bd2f4+b2ef4+c2ef4+de2f4+e3f4\n\
  cde2f4+d2e2f4+ae3f4+be3f4+ce3f4+e4f4+bc2f5+b2df5+c2df5+ad2f5+cd2f5+d3f5+acef5+adef5+bdef5+d2ef5+ce2f5+e3f5\n\
  d2e3f4+de4f4+bc2df5+abd2f5+bcd2f5+c2d2f5+bc2ef5+abdef5+bd2ef5+d3ef5+b2e2f5+ace2f5+bce2f5+cde2f5+ae3f5+de3f5\n\
";

TEST(PolyHashTable,test1) {
  std::unique_ptr<PolyRing> R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  PolyHashTable H(R.get(),3);
  std::unique_ptr<Poly> f1 = polyParseFromString(R.get(), "3bd2+7cd2+5c2f+2adf+bdf+10cef");
  PolyHashTable::MonomialArray M1, M2;
  H.fromPoly(*f1, M1);
  H.fromPoly(*f1, M2);
  EXPECT_TRUE(M2.empty());
  Poly g(*R);
  H.toPoly(M1,g);
  //  f1->display(std::cout);
  //  std::cout << std::endl;
  //  g.display(std::cout);
  //  std::cout << std::endl;
  f1->multByCoefficient(2);
  EXPECT_TRUE(g == *f1);
  //  H.dump();
  H.resize(6);
  //  H.dump();
  M1.clear();
  H.fromPoly(*f1, M1);
  Poly g2(*R);
  H.toPoly(M1,g2);
  EXPECT_TRUE(g == g2);
}

TEST(PolyHashTable,test2) {
  std::unique_ptr<PolyRing> R(ringFromString("32003 6 1\n1 1 1 1 1 1"));
  PolyHashTable H(R.get(), 3);
  std::unique_ptr<Poly> f1(polyParseFromString(R.get(), "3bd2+7cd2+5c2f+2adf+bdf+10cef"));
  std::unique_ptr<Poly> f2(polyParseFromString(R.get(), "-3bd2+4c2f+cef+f3"));
  PolyHashTable::MonomialArray M1, M2;
  H.fromPoly(*f1, M1);
  H.fromPoly(*f2, M2);
  //  H.dump(1);
}

TEST(MonomialHashTable,test1) {
  std::unique_ptr<PolyRing> R = ringFromString("32003 6 1\n1 1 1 1 1 1");
  MonomialHashTable H(R.get(), 3);
  std::unique_ptr<Poly> f1 = polyParseFromString(R.get(), "3bd2+7cd2+5c2f+2adf+bdf+10cef");
  int count = 0;
  int was_there_count = 0;
  for (int j = 0; j<10; j++)
    for (Poly::iterator i = f1->begin(); i != f1->end(); ++i)
      {
	bool was_there = H.lookupAndInsert(i.getMonomial(), count);
	count++;
	if (was_there) was_there_count++;
      }
  MonomialHashTable::Stats stats;
  H.getStats(stats);
  //H.dump(1);
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
