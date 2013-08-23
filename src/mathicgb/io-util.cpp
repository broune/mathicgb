// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "io-util.hpp"

#include "Poly.hpp"
#include "MTArray.hpp"
#include "io-util.hpp"
#include "Scanner.hpp"
#include "MathicIO.hpp"
#include "PolyHeap.hpp"
#include "PolyGeoBucket.hpp"
#include "SigPolyBasis.hpp"
#include "SignatureGB.hpp"
#include "MTArray.hpp"
#include "Basis.hpp"
#include "PolyBasis.hpp"
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>

MATHICGB_NAMESPACE_BEGIN

std::unique_ptr<Poly> polyParseFromString(const PolyRing *R, const std::string &s)
{
  std::unique_ptr<Poly> f(new Poly(*R));
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

std::unique_ptr<Basis> basisParseFromString(std::string str)
{
  std::istringstream inStream(str);
  Scanner in(inStream);
  auto p = MathicIO<>().readRing(true, in);
  auto& ring = *p.first.release(); // todo: fix leak
  return make_unique<Basis>(MathicIO<>().readBasis(ring, false, in));
}

std::unique_ptr<PolyRing> ringFromString(std::string ringinfo)
{
  std::stringstream ifil(ringinfo);
  return std::unique_ptr<PolyRing>(PolyRing::read(ifil).first);
}

Monomial stringToMonomial(const PolyRing *R, std::string mon)
{
  Monomial result = R->allocMonomial();
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

std::string toString(SigPolyBasis *I)
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

std::string toString(SigPolyBasis *I, int)
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

std::string toString(Basis *I)
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

void output(FILE* file, const PolyBasis &I)
{
  for (size_t i = 0; i < I.size(); i++)
    {
      if (!I.retired(i))
        {
          I.poly(i).display(file, false);
          fputc('\n', file);
        }
    }
}

MATHICGB_NAMESPACE_END
