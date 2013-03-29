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

TEST(MonoMonoid, ParsePrintM2) {
  MonoMonoid<int32> m(100);
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

TEST(MonoMonoid, Product) {
  typedef MonoMonoid<int32> Monoid;
  Monoid m(49);
  Monoid::MonoPool pool(m);
  auto mono = pool.alloc();
  auto checkProd = [&](const char* str) {
    auto v = parseVector(m, str);
    MATHICGB_ASSERT(v.size() == 3);
    const auto& a = v.front();
    const auto& b = *++v.begin();
    const auto& c = v.back();
    ASSERT_EQ(m.hashOfProduct(a, b), m.hash(c));
    ASSERT_EQ(m.hashOfProduct(a, b), m.hashOfProduct(b, a));

    m.multiply(a, b, mono);
    ASSERT_TRUE(m.equal(c, mono));
    ASSERT_TRUE(m.compare(c, mono) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(c), m.hash(mono));

    m.multiply(b, a, mono);
    ASSERT_TRUE(m.equal(c, mono));
    ASSERT_TRUE(m.compare(c, mono) == Monoid::EqualTo);
    ASSERT_EQ(m.hash(c), m.hash(mono));

    if (!m.isIdentity(a))
      ASSERT_TRUE(m.lessThan(b, mono));
    else
      ASSERT_TRUE(m.compare(b, mono) == Monoid::EqualTo);

    if (!m.isIdentity(b))
      ASSERT_TRUE(m.lessThan(a, mono));
    else
      ASSERT_TRUE(m.compare(a, mono) == Monoid::EqualTo);
  };
  checkProd("1 1 1");
  checkProd("a 1 a");
  checkProd("1 Vx Vx");
  checkProd("aV bx abxV");
  checkProd("a a2 a3");
  checkProd("V V2 V3");
  checkProd("arlgh svug arlg2hsvu");
  checkProd("abcdefghiV ab2c3d4e5f6g7h8i9V11 a2b3c4d5e6f7g8h9i10V12");
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
      //m.printM2(*lesser, std::cout); std::cout<<' ';
      //m.printM2(*greater, std::cout); std::cout<<std::endl;
      
      ASSERT_FALSE(m.equal(*lesser, *greater));
      ASSERT_TRUE(m.lessThan(*lesser, *greater));
      ASSERT_FALSE(m.lessThan(*greater, *lesser));
      ASSERT_EQ(m.compare(*lesser, *greater), Monoid::LessThan);
      ASSERT_EQ(m.compare(*greater, *lesser), Monoid::GreaterThan);
    }
  }
}
