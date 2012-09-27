// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include "Poly.hpp"
#include "MTArray.hpp"
#include "io-util.hpp"

#include "PolyHeap.hpp"
#include "PolyGeoBucket.hpp"
#include "GroebnerBasis.hpp"
#include "SignatureGB.hpp"
#include "MTArray.hpp"

#include "Ideal.hpp"
#include "PolyBasis.hpp"

std::unique_ptr<Poly> polyParseFromString(const PolyRing *R, const std::string &s)
{
  std::unique_ptr<Poly> f(new Poly(R));
  std::istringstream in(s);
  f->parse(in);
  return f;
}

std::string toString(const Poly *g)
{
  std::ostringstream o;
  g->display(o);
  return o.str();
}

std::unique_ptr<Ideal> idealParseFromString(std::string str)
{
  std::istringstream i(str);
  return std::unique_ptr<Ideal>(Ideal::parse(i));
}

std::unique_ptr<PolyRing> ringFromString(std::string ringinfo)
{
  std::stringstream ifil(ringinfo);
  return std::unique_ptr<PolyRing>(PolyRing::read(ifil));
}

Monomial stringToMonomial(const PolyRing *R, std::string mon)
{
  Monomial result = R->allocMonomial1();
  std::stringstream ifil(mon);
  R->monomialParse(ifil, result);
  return result;
}

std::string monomialToString(const PolyRing *R, const Monomial& mon)
{
  std::ostringstream o;
  R->monomialDisplay(o,mon);
  return o.str();
}

///////// Which of the following do we need/use? /////////////
#ifndef NEWMONOMIALS
monomial monomialFromString(const PolyRing *R, std::string mon)
{
  // This is poor code, to only be used for testing!
  monomial result = new int[R->maxMonomialSize()];
  std::stringstream ifil(mon);
  R->monomialRead(ifil, result);
  return result;
}

std::string monomialToString(const PolyRing *R, const_monomial mon)
{
  std::ostringstream o;
  R->monomialWrite(o,mon);
  return o.str();
}
#endif

monomial monomialParseFromString(const PolyRing *R, std::string mon)
{
  // This is poor code, to only be used for testing!
  monomial result = R->allocMonomial();
  std::stringstream ifil(mon);
  R->monomialParse(ifil, result);
  return result;
}

std::string monomialDisplay(const PolyRing *R, const_monomial mon)
{
  std::ostringstream o;
  R->monomialDisplay(o,mon);
  return o.str();
}
////////////////////////////////////////////////////////////////

std::string toString(GroebnerBasis *I)
{
  std::ostringstream o;
  for (size_t i=0; i<I->size(); i++)
    {
      o << "  ";
      I->poly(i).display(o, false);
      o << std::endl;
    }
  return o.str();
}

std::string toString(GroebnerBasis *I, int)
{
  std::ostringstream o;
  I->display(o);
  return o.str();
}

std::string toString(MonomialTableArray* H)
{
  std::ostringstream o;
  H->display(o, 1);
  return o.str();
}

std::string toString(Ideal *I)
{
  std::ostringstream o;
  for (size_t i=0; i<I->size(); i++)
    {
      o << "  ";
      I->getPoly(i)->display(o,false);
      o << std::endl;
    }
  return o.str();
}

void output(std::ostream &o, const PolyBasis &I)
{
  for (size_t i = 0; i < I.size(); i++)
    {
      if (!I.retired(i))
        {
          I.poly(i).display(o, false);
          o << std::endl;
        }
    }
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
