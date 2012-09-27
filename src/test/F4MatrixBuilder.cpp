#include "mathicgb/stdinc.h"

#include "mathicgb/Poly.hpp"
#include "mathicgb/PolyRing.hpp"
#include "mathicgb/F4MatrixBuilder.hpp"
#include "mathicgb/FreeModuleOrder.hpp"
#include "mathicgb/Ideal.hpp"
#include "mathicgb/PolyBasis.hpp"
#include "mathicgb/io-util.hpp"

#include <gtest/gtest.h>
#include <memory>

namespace {
  // We need a struct to keep the ring and so on alive after
  // construction - we cannot just return an object.
  //
  // @todo: This whole thing is fairly ridiculous - some kind of more
  // general dependency injection mechanism might be nice here.
  struct BuilderMaker {
    BuilderMaker():
      mRing(ringFromString("101 6 1\n1 1 1 1 1 1")),
      mIdeal(*mRing),
      mOrder(FreeModuleOrder::makeOrder(1, &mIdeal)),
      mBasis(*mRing, *mOrder, DivisorLookup::makeFactory(*mRing, 1)->create(true, true)) {
    }

    const Poly& addBasisElement(const std::string& str) {
      std::unique_ptr<Poly> p(new Poly(mRing.get()));
      std::istringstream in(str);
      p->parse(in);
      mBasis.insert(std::move(p));
      return mBasis.poly(mBasis.size() - 1);
    }

    F4MatrixBuilder& create() {
      MATHICGB_ASSERT(mBuilder.get() == 0);
      mBuilder.reset(new F4MatrixBuilder(mBasis));
      return *mBuilder;
    }

    const PolyRing& ring() const {return *mRing;}
     
  private:
    std::unique_ptr<PolyRing> mRing;
    Ideal mIdeal;
    std::unique_ptr<FreeModuleOrder> mOrder;
    PolyBasis mBasis;
    std::unique_ptr<F4MatrixBuilder> mBuilder;
  };
}

TEST(F4MatrixBuilder, Empty) {
  BuilderMaker maker;
  F4MatrixBuilder& builder = maker.create();

  QuadMatrix matrix;
  builder.buildMatrixAndClear(matrix);
  ASSERT_EQ(0, matrix.topLeft.rowCount());
  ASSERT_EQ(0, matrix.bottomLeft.rowCount());
  ASSERT_EQ(0, matrix.topLeft.colCount());
  ASSERT_EQ(0, matrix.topRight.colCount());
  ASSERT_EQ(0, matrix.leftColumnMonomials.size());
  ASSERT_EQ(0, matrix.rightColumnMonomials.size());
}

TEST(F4MatrixBuilder, OneByOne) {
  BuilderMaker maker;
  const Poly& p = maker.addBasisElement("a");
  F4MatrixBuilder& builder = maker.create();
  builder.addRowToMatrix(p.getLeadMonomial(), p);
  QuadMatrix qm;
  builder.buildMatrixAndClear(qm);
  const char* str = 
    "Left columns: a2\n"
    "Right columns:\n"
    "0: 0#1              | 0:                 \n"
    "                    |                    \n"
    "matrix with no rows | matrix with no rows\n";
  ASSERT_EQ(str, qm.toString());
}

TEST(F4MatrixBuilder, DirectReducers) {
  BuilderMaker maker;
  maker.addBasisElement("a6<0>"); // reducer ==, but won't be used as not added
  maker.addBasisElement("a3b2<0>+a3c"); // reducer ==
  maker.addBasisElement("c<0>"); // reducer divides
  maker.addBasisElement("d2<0>"); // does not divide
  F4MatrixBuilder& builder = maker.create();

  Poly p1(&builder.ring());
  { 
    std::istringstream in("a3<0>+b2+c+d");
    p1.parse(in);
    builder.addRowToMatrix(p1.getLeadMonomial(), p1);
  }

  Poly p2(&builder.ring());
  {
    std::istringstream in("a3<0>+2b2+3c+4d");
    p2.parse(in);
    builder.addRowToMatrix(p2.getLeadMonomial(), p2);
  }

  QuadMatrix qm;
  builder.buildMatrixAndClear(qm);

  const char* str =
    "Left columns: a6 a3b2 a3c\n"
    "Right columns: a3d\n"
    "0: 2#1         | 0:    \n"
    "1: 1#1 2#1     | 1:    \n"
    "2: 0#1 1#1 2#1 | 2: 0#1\n"
    "               |       \n"
    "0: 0#1 1#2 2#3 | 0: 0#4\n";
  // This quest is currently fragile because of the possibility of
  // reordering of the rows.
  ASSERT_EQ(str, qm.toString());
}

TEST(F4MatrixBuilder, IteratedReducer) {
  BuilderMaker maker;
  const Poly& p1 = maker.addBasisElement("a4-a3");
  const Poly& p2 = maker.addBasisElement("a-1");
  F4MatrixBuilder& builder = maker.create();
  builder.addRowToMatrix(p1.getLeadMonomial(), p2);
  QuadMatrix qm;
  builder.buildMatrixAndClear(qm);
  const char* str = 
    "Left columns: a5 a4 a3 a2 a\n"
    "Right columns: 1\n"
    "0: 0#1 1#100        | 0:                 \n"
    "1: 1#1 2#100        | 1:                 \n"
    "2: 2#1 3#100        | 2:                 \n"
    "3: 3#1 4#100        | 3:                 \n"
    "4: 4#1              | 4: 0#100           \n"
    "                    |                    \n"
    "matrix with no rows | matrix with no rows\n";
  ASSERT_EQ(str, qm.toString());
}
