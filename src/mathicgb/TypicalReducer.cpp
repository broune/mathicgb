#include "stdinc.h"
#include "TypicalReducer.hpp"

#include "SigPolyBasis.hpp"
#include "PolyBasis.hpp"
#include <iostream>

size_t TypicalReducer::preferredSetSize() const {
  return 1;
}

void TypicalReducer::reset()
{
  mArena.freeAllAllocs();
  resetReducer();
}

size_t TypicalReducer::getMemoryUse() const {
  return mArena.getMemoryUse();
}

Poly* TypicalReducer::regularReduce(
  const_monomial sig,
  const_monomial multiple,
  size_t basisElement,
  const SigPolyBasis& basis)
{
  const PolyRing& ring = basis.ring();
  ++mSigStats.reductions;

  monomial tproduct = ring.allocMonomial(mArena);
  monomial u = ring.allocMonomial(mArena);
  ring.monomialMult(multiple, basis.getLeadMonomial(basisElement), tproduct);

  size_t reducer = basis.regularReducer(sig, tproduct);
  if (reducer == static_cast<size_t>(-1)) {
    ++mSigStats.singularReductions;
    mArena.freeAllAllocs();
    return 0; // singular reduction: no regular top reduction possible
  }

  ring.monomialDivide(tproduct, basis.getLeadMonomial(reducer), u);

  coefficient coef;
  ring.coefficientSet(coef, 1);
  insertTail(const_term(coef, multiple), &basis.poly(basisElement));

  MATHICGB_ASSERT(ring.coefficientIsOne(basis.getLeadCoefficient(reducer)));
  ring.coefficientFromInt(coef, -1);
  insertTail(const_term(coef, u), &basis.poly(reducer));
  basis.basis().usedAsReducer(reducer);

  Poly* result = new Poly(ring);

  unsigned long long steps = 2; // number of steps in this reduction
  for (const_term v; leadTerm(v);) {
    MATHICGB_ASSERT(v.coeff != 0);
    reducer = basis.regularReducer(sig, v.monom);
    if (reducer == static_cast<size_t>(-1)) { // no reducer found
      result->appendTerm(v.coeff, v.monom);
      removeLeadTerm();
    } else { // reduce by reducer
      ++steps;
      basis.basis().usedAsReducer(reducer);
      monomial mon = ring.allocMonomial(mArena);
      ring.monomialDivide(v.monom, basis.getLeadMonomial(reducer), mon);
      ring.coefficientDivide(v.coeff, basis.getLeadCoefficient(reducer), coef);
      ring.coefficientNegateTo(coef);
      removeLeadTerm();
      insertTail(const_term(coef, mon), &basis.poly(reducer));
    }
  }
  result->makeMonic();

  mSigStats.steps += steps;
  mSigStats.maxSteps = std::max(mSigStats.maxSteps, steps);
  if (result->isZero())
    ++mSigStats.zeroReductions;

  reset();
  return result;
}

std::unique_ptr<Poly> TypicalReducer::classicReduce(const Poly& poly, const PolyBasis& basis) {
  monomial identity = basis.ring().allocMonomial(mArena);
  basis.ring().monomialSetIdentity(identity);
  insert(identity, &poly);

  return classicReduce(basis);
}

std::unique_ptr<Poly> TypicalReducer::classicTailReduce(const Poly& poly, const PolyBasis& basis) {
  MATHICGB_ASSERT(&poly.ring() == &basis.ring());
  MATHICGB_ASSERT(!poly.isZero());
  term identity;
  identity.monom = basis.ring().allocMonomial(mArena);
  basis.ring().monomialSetIdentity(identity.monom);
  basis.ring().coefficientSetOne(identity.coeff);
  insertTail(identity, &poly);

  std::unique_ptr<Poly> result(new Poly(basis.ring()));
  result->appendTerm(poly.getLeadCoefficient(), poly.getLeadMonomial());

  return classicReduce(std::move(result), basis);
}

std::unique_ptr<Poly> TypicalReducer::classicReduceSPoly(
  const Poly& a,
  const Poly& b,
  const PolyBasis& basis
) {
  const PolyRing& ring = basis.ring();

  monomial lcm = ring.allocMonomial();
  ring.monomialLeastCommonMultiple
    (a.getLeadMonomial(), b.getLeadMonomial(), lcm);

  // insert tail of multiple of a
  monomial multiple1 = ring.allocMonomial();
  ring.monomialDivide(lcm, a.getLeadMonomial(), multiple1);
  coefficient plusOne;
  ring.coefficientSet(plusOne, 1);
  insertTail(const_term(plusOne, multiple1), &a);

  // insert tail of multiple of b
  monomial multiple2 = ring.allocMonomial();
  ring.monomialDivide(lcm, b.getLeadMonomial(), multiple2);
  coefficient minusOne = plusOne;
  ring.coefficientNegateTo(minusOne);
  insertTail(const_term(minusOne, multiple2), &b);

  std::unique_ptr<Poly> reduced = classicReduce(basis);
  ring.freeMonomial(lcm);
  ring.freeMonomial(multiple1);
  ring.freeMonomial(multiple2);
  return std::move(reduced);
}

void TypicalReducer::classicReduceSPolySet
(std::vector<std::pair<size_t, size_t> >& spairs,
 const PolyBasis& basis,
 std::vector<std::unique_ptr<Poly> >& reducedOut) {
  for (auto it = spairs.begin(); it != spairs.end(); ++it) {
    auto reducedSPoly =
      classicReduceSPoly(basis.poly(it->first), basis.poly(it->second), basis);
    if (!reducedSPoly->isZero())
      reducedOut.push_back(std::move(reducedSPoly));
  }
}

void TypicalReducer::classicReducePolySet
(const std::vector<std::unique_ptr<Poly> >& polys,
 const PolyBasis& basis,
 std::vector<std::unique_ptr<Poly> >& reducedOut)
{
  for (auto it = polys.begin(); it != polys.end(); ++it) {
    auto reducedPoly = classicReduce(**it, basis);
    if (!reducedPoly->isZero())
      reducedOut.push_back(std::move(reducedPoly));
  }  
}

void TypicalReducer::setMemoryQuantum(size_t quantum) {
}

std::unique_ptr<Poly> TypicalReducer::classicReduce
    (std::unique_ptr<Poly> result, const PolyBasis& basis) {
  const PolyRing& ring = basis.ring();
  MATHICGB_ASSERT(&result->ring() == &ring);
  ++mClassicStats.reductions;

  if (tracingLevel > 100)
    std::cerr << "Classic reduction begun." << std::endl;

  coefficient coef;
  unsigned long long steps = 0; // number of steps in this reduction
  for (const_term v; leadTerm(v);) {
    if (tracingLevel > 100) {
      std::cerr << "from reducer queue: ";
      basis.ring().monomialDisplay(std::cerr, v.monom);
      std::cerr << std::endl;
    }

    size_t reducer = basis.classicReducer(v.monom);
    if (reducer == static_cast<size_t>(-1)) { // no reducer found
      MATHICGB_ASSERT(
        result->isZero() ||
        basis.monoid().lessThan(v.monom, result->backMonomial())
      );
      result->appendTerm(v.coeff, v.monom);
      removeLeadTerm();
    } else { // reduce by reducer
      ++steps;
      basis.usedAsReducer(reducer);
      monomial mon = ring.allocMonomial(mArena);
      ring.monomialDivide(v.monom, basis.leadMonomial(reducer), mon);
      ring.coefficientDivide(v.coeff, basis.leadCoefficient(reducer), coef);
      ring.coefficientNegateTo(coef);
      removeLeadTerm();
      insertTail(const_term(coef, mon), &basis.poly(reducer));

      if (tracingLevel > 100) {
        std::cerr << "Reducing by basis element " << reducer << ": ";
        basis.poly(reducer).display(std::cerr);
        std::cerr << std::endl;
        std::cerr << "multiplied by: " << coef << "  *  ";
        basis.ring().monomialDisplay(std::cerr, mon);
        std::cerr << std::endl;
      }
    }
  }
  result->makeMonic();

  mClassicStats.steps += steps;
  mClassicStats.maxSteps = std::max(mClassicStats.maxSteps, steps);
  if (result->isZero())
    ++mClassicStats.zeroReductions;

  if (tracingLevel > 100)
    std::cerr << "Classic reduction done." << std::endl;

  reset();
  return std::move(result);
}

std::unique_ptr<Poly> TypicalReducer::classicReduce(const PolyBasis& basis) {
  return classicReduce(make_unique<Poly>(basis.ring()), basis);
}
