#include "mathicgb/stdinc.h"

#include "mathicgb/MonoMonoid.hpp"
#include <gtest/gtest.h>
#include <sstream>

// expect(i,j) encodes a matrix with interesting bit patterns that
// are supposed to be likely to surface errors in how monomials are
// stored inside a vector.
uint32 expect(size_t mono, size_t var, size_t varCount) {
  const auto unique = var + varCount * mono + 1;

  while (true) {
    // 000
    if (mono == 0)
      return 0;
    --mono;

    // 100
    // 010
    // 001
    if (mono < varCount)
      return var == mono ? unique : 0;
    mono -= varCount;

    // 000
    // 100
    // 110
    // 111
    if (mono < varCount + 1)
      return var < mono ? unique : 0;
    mono -= varCount + 1;

    // 111
    // 011
    // 001
    // 000
    if (mono < varCount + 1)
      return var >= mono ? unique : 0;
    mono -= varCount + 1;

    // 101
    // 010
    if (mono < 4)
      return (var % 2) == (mono % 2) ? unique : 0;
    mono -= 4;

    // 100100
    // 010010
    // 001001
    if (mono < 6)
      return (var % 3) == (mono % 3) ? unique : 0;
    mono -= 6;

    // mix the patterns
    mono += var % 17;
  }
};

TEST(MonoMonoid, VarCount) {
  ASSERT_EQ(0, MonoMonoid<int8>(0).varCount());
  ASSERT_EQ(1000 * 1000, MonoMonoid<int8>(1000 * 1000).varCount());
  ASSERT_EQ(1, MonoMonoid<int16>(1).varCount());
  ASSERT_EQ(2, MonoMonoid<int32>(2).varCount());
  ASSERT_EQ(12, MonoMonoid<int64>(12).varCount());
}

TEST(MonoMonoid, MonoVector) {
  typedef MonoMonoid<int32> Monoid;
  typedef Monoid::VarIndex VarIndex;
  typedef Monoid::MonoVector MonoVector;

  Monoid monoid(13);
  MonoVector v(monoid);
  MonoVector v2(monoid);
  ASSERT_EQ(v2.monoid(), monoid);
  const auto varCount = monoid.varCount();

  ASSERT_TRUE(v.empty());
  size_t count = 1000;


  // Not a correctness error, but empty vectors should preferably not
  // use any memory.
  ASSERT_EQ(0, v.memoryBytesUsed());

  for (size_t i = 0; i < count; ++i) {
    ASSERT_EQ(i, v.size());
    v.push_back(); // push_back, no param
    ASSERT_GT(v.memoryBytesUsed(), 0);
    ASSERT_FALSE(v.empty()); // empty
    ASSERT_EQ(i + 1, v.size()); // size


    ASSERT_TRUE(monoid.isIdentity(v.back())); // isIdentity true, back non-const
    bool allZero = true;
    for (VarIndex var = 0; var < varCount; ++var) {
      const auto exponent = expect(i, var, varCount);
      if (exponent != 0) {
	allZero = false;
	monoid.setExponent(var, exponent, v.back());
      }
    }
    ASSERT_EQ(allZero, monoid.isIdentity(v.back())); // isIdentity false
    v2.push_back(v.back()); // push_back with param
    ASSERT_TRUE(monoid.equal(v.back(), v2.back()));
  }
  auto it = v.begin();
  ASSERT_EQ(it, v.cbegin());
  for (size_t i = 0; i < count; ++i, ++it) {
    MonoVector::const_iterator tmp;
    ASSERT_TRUE(it != tmp);
    tmp = it;
    ASSERT_EQ(tmp, it);
    ASSERT_TRUE(v.end() != it);

    for (VarIndex var = 0; var < monoid.varCount(); ++var) {
      ASSERT_EQ(expect(i, var, varCount), monoid.exponent(*it, var));
    }
  }
  ASSERT_EQ(v.end(), it);
  ASSERT_EQ(v.cend(), it);

  ASSERT_EQ(v, v2); // operator== true
  monoid.setExponent(0, 1 + monoid.exponent(v2.back(), 0), v2.back());
  ASSERT_TRUE(v != v2); // operator!=, true, same length

  auto& vc = const_cast<const MonoVector&>(v);
  ASSERT_TRUE(monoid.equal(v.front(), *v2.begin())); // front, non-const
  ASSERT_TRUE(monoid.equal(vc.front(), *v2.begin())); // front, const
  ASSERT_TRUE(monoid.equal(vc.back(), v.back())); // back, non-const

  auto v3(v2); // copy constructor
  ASSERT_EQ(v3.monoid(), monoid);
  ASSERT_TRUE(v != v3 && v2 == v3);
  v2.swap(v); // member swap
  ASSERT_TRUE(v == v3 && v2 != v3);
  std::swap(v, v2); // std::swap
  ASSERT_TRUE(v != v3 && v2 == v3);
  using std::swap;
  swap(v, v2); // let compiler decide which swap to use
  ASSERT_TRUE(v == v3 && v2 != v3);
  swap(v, v2); // get back to original state
  ASSERT_TRUE(v != v3 && v2 == v3);

  ASSERT_FALSE(v3 != v2); // operator!=, false, same length
  v3.push_back();
  ASSERT_TRUE(v3 != v2); // operator!=, true, different length
  

  ASSERT_FALSE(v3 == v);
  v3 = v; // copy assignment
  ASSERT_EQ(v3.monoid(), monoid);
  ASSERT_EQ(v3, v);

  ASSERT_FALSE(v3.empty());
  v2 = std::move(v3); // move assignment
  ASSERT_EQ(v2.monoid(), monoid);
  ASSERT_EQ(v2, v);
  ASSERT_TRUE(v3.empty());

  ASSERT_FALSE(v2.empty());
  auto v4(std::move(v2)); // move constructor
  ASSERT_EQ(v4.monoid(), monoid);
  ASSERT_TRUE(v2.empty());
  ASSERT_EQ(v4, v);

  ASSERT_FALSE(v.empty());
  v.clear();
  ASSERT_TRUE(v.empty());
}

TEST(MonoMonoid, MonoPool) {
  typedef MonoMonoid<int32> Monoid;
  typedef Monoid::VarIndex VarIndex;
  typedef Monoid::Mono Mono;

  for (int q = 0; q < 2; ++q) {
    Monoid monoid(13);
    Monoid::MonoPool pool(monoid);
    const auto varCount = monoid.varCount();

    const auto count = 1000;
    std::vector<Mono> monos;
    for (int i = 0; i < count; ++i) {
      pool.alloc();
      pool.free(pool.alloc());
      auto m1 = pool.alloc();
      ASSERT_TRUE(monoid.isIdentity(m1));
      auto m2 = pool.alloc();
      ASSERT_TRUE(monoid.isIdentity(m2));
      for (VarIndex var = 0; var < varCount; ++var) {
	monoid.setExponent(var, 1, m1);
	monoid.setExponent(var, 1, m2);
      }
      if (i > 10) {
	using std::swap;
	swap(m2, monos[i - 10]);
      }
      monos.push_back(std::move(m1));
    }

    // This ensures that we get to each entry in monos exactly once.
    MATHICGB_ASSERT((count % 17) != 0); 
    int i = 0;
    do {
      MATHICGB_ASSERT(!monos[i].isNull());
      ASSERT_FALSE(monoid.isIdentity(monos[i]));
      pool.free(monos[i]);
      ASSERT_TRUE(monos[i].isNull());
      pool.free(monos[i]);
      ASSERT_TRUE(monos[i].isNull());
      i = (i + 17) % count;
    } while (i != 0);

    // If the ordering of monomials inside the pool has anything to do with
    // allocation and deallocation order, then the monomials inside the
    // pool are at this point all jumbled around. All the entries were also
    // non-zero before, so we test that new allocations are the identity.

    for (int i = 0; i < count; ++i) {
      monos[i] = pool.alloc();
      ASSERT_TRUE(monoid.isIdentity(monos[i]));
      for (VarIndex var = 0; var < varCount; ++var)
	monoid.setExponent(var, expect(i, var, varCount), monos[i]);
    }
    for (int i = 0; i < count; ++i) {
      for (VarIndex var = 0; var < varCount; ++var) {
	ASSERT_EQ(expect(i, var, varCount), monoid.exponent(monos[i], var));
      }
    }
    // everything should be free'd now. Let's do all that again.
  }
}

namespace {
  template<class M>
  typename M::MonoVector parseVector(M& monoid, const char* str) {
    typename M::MonoVector v(monoid);
    std::istringstream in(str);
    v.parseM2(in);
    return v;
  }
}

namespace {
  template<class E>
  void parsePrintM2Helper() {
    MonoMonoid<E> m(100);
    const char* str = "1 a z A Z ab a2 a2b ab2 a20b30 1<1> a<2> a2<3> ab<11>\n";
    auto v2 = parseVector(m, str);
    std::ostringstream v2Out;
    v2.printM2(v2Out);
    ASSERT_EQ(str, v2Out.str());

    decltype(v2) v(m);
    v.push_back(); // 1

    v.push_back(); // a
    m.setExponent(0, 1, v.back());
 
    v.push_back(); // z
    m.setExponent(25, 1, v.back());

    v.push_back(); // A
    m.setExponent(26, 1, v.back());

    v.push_back(); // Z
    m.setExponent(51, 1, v.back());

    v.push_back(); // ab
    m.setExponent(0, 1, v.back());
    m.setExponent(1, 1, v.back());

    v.push_back(); // a2
    m.setExponent(0, 2, v.back());

    v.push_back(); // a2b
    m.setExponent(0, 2, v.back());
    m.setExponent(1, 1, v.back());

    v.push_back(); // ab2
    m.setExponent(0, 1, v.back());
    m.setExponent(1, 2, v.back());

    v.push_back(); // a20b30
    m.setExponent(0, 20, v.back());
    m.setExponent(1, 30, v.back());

    v.push_back(); // 1<1>
    m.setComponent(1, v.back());

    v.push_back(); // a<2>
    m.setComponent(2, v.back());
    m.setExponent(0, 1, v.back());

    v.push_back(); // a2<3>
    m.setComponent(3, v.back());
    m.setExponent(0, 2, v.back());

    v.push_back(); // ab<11>
    m.setComponent(11, v.back());
    m.setExponent(0, 1, v.back());
    m.setExponent(1, 1, v.back());

    std::ostringstream vOut;
    v.printM2(vOut);
    ASSERT_EQ(str, vOut.str());
  
    ASSERT_EQ(v, v2);
  }
}

TEST(MonoMonoid, ParsePrintM2) {
  parsePrintM2Helper<int32>();
  parsePrintM2Helper<int16>();
  parsePrintM2Helper<int8>();
}

TEST(MonoMonoid, MultiplyDivide) {
  typedef MonoMonoid<int32> Monoid;
  Monoid m(49);
  Monoid::MonoPool pool(m);
  auto mono = pool.alloc();
  auto check = [&](const char* str) {
    auto v = parseVector(m, str);
    MATHICGB_ASSERT(v.size() == 3);
    const auto& a = v.front();
    const auto& b = *++v.begin();
    const auto& c = v.back();
    ASSERT_EQ(m.hashOfProduct(a, b), m.hash(c));
    ASSERT_EQ(m.hashOfProduct(a, b), m.hashOfProduct(b, a));

    // isProductOf
    ASSERT_TRUE(m.isProductOf(a, b, c));
    ASSERT_TRUE(m.isProductOfHintTrue(a, b, c));
    ASSERT_TRUE(m.isTwoProductsOfHintTrue(a, a, b, c, c));
    

    // a*b == c using multiply
    m.multiply(a, b, mono);
    ASSERT_TRUE(m.equal(c, mono));
    ASSERT_TRUE(m.compare(c, mono) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(c), m.hash(mono));

    // c/a == b using divide
    m.divide(a, c, mono);
    ASSERT_TRUE(m.equal(b, mono));
    ASSERT_TRUE(m.compare(b, mono) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(b), m.hash(mono));

    // c/b == a using divideInPlace
    m.copy(c, mono);
    m.divideInPlace(b, mono);
    ASSERT_TRUE(m.equal(a, mono));
    ASSERT_TRUE(m.compare(a, mono) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(a), m.hash(mono));

    // a*b == c using multiplyInPlace
    m.copy(a, mono);
    m.multiplyInPlace(b, mono);
    ASSERT_TRUE(m.equal(c, mono));
    ASSERT_TRUE(m.compare(c, mono) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(c), m.hash(mono));

    // check properties that mono=a*b should have
    ASSERT_TRUE(m.divides(mono, c));
    ASSERT_TRUE(m.divides(c, mono));
    ASSERT_TRUE(m.divides(a, mono));
    ASSERT_TRUE(m.divides(b, mono));

    if (!m.isIdentity(a)) {
      ASSERT_TRUE(m.lessThan(b, mono));
      ASSERT_FALSE(m.lessThan(mono, b));
      ASSERT_TRUE(m.compare(mono, b) == Monoid::GreaterThan);
      ASSERT_FALSE(m.divides(mono, b));

      ASSERT_FALSE(m.isProductOf(a, c, b));
      ASSERT_FALSE(m.isProductOfHintTrue(a, c, b));
      ASSERT_FALSE(m.isTwoProductsOfHintTrue(c, c, a, b, b));
      ASSERT_FALSE(m.isTwoProductsOfHintTrue(b, c, a, c, b));
      ASSERT_FALSE(m.isTwoProductsOfHintTrue(c, b, a, b, c));
    } else {
      ASSERT_TRUE(m.equal(b, mono));
      ASSERT_TRUE(m.compare(b, mono) == Monoid::EqualTo);
      ASSERT_TRUE(m.divides(mono, b));
    }

    if (!m.isIdentity(b)) {
      ASSERT_TRUE(m.lessThan(a, mono));
      ASSERT_FALSE(m.lessThan(mono, a));
      ASSERT_TRUE(m.compare(mono, a) == Monoid::GreaterThan);
      ASSERT_FALSE(m.divides(mono, a));

      ASSERT_FALSE(m.isProductOf(c, b, a));
      ASSERT_FALSE(m.isProductOfHintTrue(b, c, a));
      ASSERT_FALSE(m.isTwoProductsOfHintTrue(c, c, b, a, a));
      ASSERT_FALSE(m.isTwoProductsOfHintTrue(a, c, b, c, a));
      ASSERT_FALSE(m.isTwoProductsOfHintTrue(c, a, b, a, c));
    } else {
      ASSERT_TRUE(m.equal(a, mono));
      ASSERT_TRUE(m.compare(a, mono) == Monoid::EqualTo);
      ASSERT_TRUE(m.divides(mono, a));
    }

    // Check that aliased parameters work.
    m.multiply(mono, mono, mono);
    m.divide(mono, mono, mono);
    MATHICGB_ASSERT(m.isIdentity(mono));

    // Check that negative exponents work.
    m.divideToNegative(a, b, mono);
    m.multiply(a, mono, mono);
    ASSERT_TRUE(m.equal(mono, b));
    
    m.divideToNegative(b, a, mono);
    m.multiply(b, mono, mono);
    ASSERT_TRUE(m.equal(mono, a));
  };
  check("1 1 1");
  check("a<5> 1 a<5>");
  check("1 Vx Vx");
  check("aV bx abxV");
  check("a a2 a3");
  check("V<2> V2 V3<2>");
  check("arlgh svug arlg2hsvu");
  check("abcdefghiV<7> ab2c3d4e5f6g7h8i9V11 a2b3c4d5e6f7g8h9i10V12<7>");
}

TEST(MonoMonoid, LcmColon) {
  typedef MonoMonoid<int32> Monoid;
  Monoid m(49);
  Monoid::MonoPool pool(m);
  auto mono = pool.alloc();
  auto mono2 = pool.alloc();
  auto check = [&](const char* str) {
    auto v = parseVector(m, str);
    MATHICGB_ASSERT(v.size() == 3);
    const auto& a = v.front();
    const auto& b = *++v.begin();
    const auto& lcm = v.back();

    // isLcm
    ASSERT_TRUE(m.isLcm(a, b, lcm));
    m.copy(lcm, mono);
    m.setExponent(1, m.exponent(mono, 1) + 1, mono);
    ASSERT_FALSE(m.isLcm(a, b, mono));

    // dividesLcm
    ASSERT_TRUE(m.dividesLcm(lcm, a, b));
    ASSERT_FALSE(m.dividesLcm(mono, a, b));
    ASSERT_TRUE(m.dividesLcm(a, a, a));
    ASSERT_TRUE(m.dividesLcm(a, a, b));
    ASSERT_TRUE(m.dividesLcm(b, b, b));
    ASSERT_TRUE(m.dividesLcm(b, b, a));

    // lcm(a, b)
    m.lcm(a, b, mono);
    ASSERT_TRUE(m.equal(mono, lcm));
    ASSERT_TRUE(m.compare(mono, lcm) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(lcm), m.hash(mono));

    // lcm(b, a)
    m.lcm(b, a, mono);
    ASSERT_TRUE(m.equal(mono, lcm));
    ASSERT_TRUE(m.compare(mono, lcm) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(lcm), m.hash(mono));

    // colons
    m.colons(a, b, mono, mono2);
    m.multiply(b, mono, mono);
    m.multiply(a, mono2, mono2);
    ASSERT_TRUE(m.equal(mono, lcm));
    ASSERT_TRUE(m.compare(mono, lcm) == Monoid::EqualTo);
    ASSERT_TRUE(m.equal(mono2, lcm));
    ASSERT_TRUE(m.compare(mono2, lcm) == Monoid::EqualTo);
  };
  check("1 1 1");
  check("a<2> 1<2> a<2>");
  check("1 Vx Vx");
  check("aV bx abxV");
  check("a a2 a2");
  check("V<3> V2<3> V2<3>");
  check("arlgh svug arlghsvu");
  check("a6b7c8d9efghiV ab2c3d4e5f6g7h8i9V11 a6b7c8d9e5f6g7h8i9V11");
}

TEST(MonoMonoid, Order) {
  typedef MonoMonoid<int32> Monoid;
  Monoid m(52);
  auto v = parseVector(m, "1 Z A z c b a c2 bc ac b2 ab a2 c3 abc b3 a3");

  for (auto greater = v.begin(); greater != v.end(); ++greater) {
    ASSERT_EQ(m.compare(*greater, *greater), Monoid::EqualTo);
    ASSERT_TRUE(m.equal(*greater, *greater));
    ASSERT_FALSE(m.lessThan(*greater, *greater));

    for (auto lesser = v.begin(); lesser != greater; ++lesser) {
      ASSERT_FALSE(m.equal(*lesser, *greater));
      ASSERT_TRUE(m.lessThan(*lesser, *greater));
      ASSERT_FALSE(m.lessThan(*greater, *lesser));
      ASSERT_EQ(m.compare(*lesser, *greater), Monoid::LessThan);
      ASSERT_EQ(m.compare(*greater, *lesser), Monoid::GreaterThan);
    }
  }
}

TEST(MonoMonoid, RelativelyPrime) {
  typedef MonoMonoid<int32> Monoid;
  Monoid m(49);
  Monoid::MonoPool pool(m);
  auto mono = pool.alloc();
  auto mono2 = pool.alloc();
  auto check = [&](const char* str, bool relativelyPrime) {
    auto v = parseVector(m, str);
    MATHICGB_ASSERT(v.size() == 2);
    ASSERT_EQ(relativelyPrime, m.relativelyPrime(v.front(), v.back()));
    ASSERT_EQ(relativelyPrime, m.relativelyPrime(v.back(), v.front()));
  };
  check("1 1", true);
  check("1 abcdefgh", true);
  check("abc defgh", true);
  check("bdfh aceg", true);
  check("bdefh aceg", false);
  check("abcdefgh abcdefgh", false);
  check("fgh abcdef", false);
}

TEST(MonoMonoid, SetExponents) {
  typedef MonoMonoid<int32> Monoid;
  Monoid m(5);
  auto v = parseVector(m, "a1b2c3d4e5");
  int32 exponents[] = {1, 2, 3, 4, 5};
  v.push_back();
  m.setExponents(exponents, v.back());
  ASSERT_TRUE(m.equal(v.front(), v.back()));  
}

TEST(MonoMonoid, HasAmpleCapacity) {
  // non-total-degree grading, last char is c
  std::vector<int8> v8;
  v8.push_back(1);
  v8.push_back(10);
  v8.push_back(1);
  MonoMonoid<int8> m8(v8);

  // total degree grading, last char is d.
  MonoMonoid<int16> m16(4);

  // non-total-degree grading, last char is e
  std::vector<int32> v32;
  v32.push_back(1);
  v32.push_back(10);
  v32.push_back(1);
  v32.push_back(1);
  v32.push_back(1);
  MonoMonoid<int32> m32(v32);

  // total degree grading, last char is n
  MonoMonoid<int32> m32t(14);


  // pure power, first variable
  auto f8 = parseVector(m8, "a63 a64");
  ASSERT_TRUE(m8.hasAmpleCapacity(f8.front()));
  ASSERT_FALSE(m8.hasAmpleCapacity(f8.back()));

  auto f16 = parseVector(m16, "a16383 a16384");
  ASSERT_TRUE(m16.hasAmpleCapacity(f16.front()));
  ASSERT_FALSE(m16.hasAmpleCapacity(f16.back()));

  auto f32 = parseVector(m32, "a1073741823 a1073741824");
  ASSERT_TRUE(m32.hasAmpleCapacity(f32.front()));
  ASSERT_FALSE(m32.hasAmpleCapacity(f32.back()));

  auto f32t = parseVector(m32t, "a1073741823 a1073741824");
  ASSERT_TRUE(m32.hasAmpleCapacity(f32t.front()));
  ASSERT_FALSE(m32.hasAmpleCapacity(f32t.back()));

  // pure power, last variable
  auto l8 = parseVector(m8, "c63 c64");
  ASSERT_TRUE(m8.hasAmpleCapacity(l8.front()));
  ASSERT_FALSE(m8.hasAmpleCapacity(l8.back()));

  auto l16 = parseVector(m16, "d16383 d16384");
  ASSERT_TRUE(m16.hasAmpleCapacity(l16.front()));
  ASSERT_FALSE(m16.hasAmpleCapacity(l16.back()));

  auto l32 = parseVector(m32, "e1073741823 e1073741824");
  ASSERT_TRUE(m32.hasAmpleCapacity(l32.front()));
  ASSERT_FALSE(m32.hasAmpleCapacity(l32.back()));

  auto l32t = parseVector(m32t, "n1073741823 n1073741824");
  ASSERT_TRUE(m32t.hasAmpleCapacity(l32t.front()));
  ASSERT_FALSE(m32t.hasAmpleCapacity(l32t.back()));

  // no exponent is too high but the degree is
  auto d8 = parseVector(m8, "abc52 abc53");
  ASSERT_TRUE(m8.hasAmpleCapacity(d8.front()));
  ASSERT_FALSE(m8.hasAmpleCapacity(d8.back()));

  auto d16 = parseVector(m16, "abcd16380 abcd16381");
  ASSERT_TRUE(m16.hasAmpleCapacity(d16.front()));
  ASSERT_FALSE(m16.hasAmpleCapacity(d16.back()));

  auto d32 = parseVector(m32, "abcde1073741810 abcde1073741811");
  ASSERT_TRUE(m32.hasAmpleCapacity(d32.front()));
  ASSERT_FALSE(m32.hasAmpleCapacity(d32.back()));

  auto d32t = parseVector
    (m32t, "abcdefghijklmn1073741810 abcdefghijklmn1073741811");
  ASSERT_TRUE(m32t.hasAmpleCapacity(d32t.front()));
  ASSERT_FALSE(m32t.hasAmpleCapacity(d32t.back()));
}
