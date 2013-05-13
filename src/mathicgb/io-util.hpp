// Copyright 2011 Michael E. Stillman

#ifndef _io_util_h_
#define _io_util_h_

#include "PolyRing.hpp"

class Poly;
class SigPolyBasis;
class MonomialTableArray;
class PolyBasis;
class Basis;

std::unique_ptr<PolyRing> ringFromString(std::string ringinfo);
monomial monomialFromString(const PolyRing *R, std::string mon);
monomial monomialParseFromString(const PolyRing *R, std::string mon);
std::string monomialToString(const PolyRing *R, const_monomial mon);
std::string monomialDisplay(const PolyRing *R, const_monomial mon);

Monomial stringToMonomial(const PolyRing *R, std::string mon);
std::string monomialToString(const PolyRing *R, const Monomial& mon);

std::string toString(SigPolyBasis *);
std::string toString(MonomialTableArray *);
std::string toString(SigPolyBasis *, int unused); // also displays signature
std::string toString(Basis *);
std::string toString(const Poly *);

std::unique_ptr<Basis> basisParseFromString(std::string str);
std::unique_ptr<Poly> polyParseFromString(const PolyRing *R, const std::string &s);

void output(std::ostream &o, const PolyBasis &I);
void output(FILE* file, const PolyBasis &I);

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
