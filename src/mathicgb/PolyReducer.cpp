// Copyright 2011 Michael E. Stillman
#include "stdinc.h"
#include "PolyReducer.hpp"

PolyReducer::PolyReducer(const PolyRing *R0):
  R(R0), mMemUsage(0)
{
  f = new Poly(*R);
  f_iter = f->begin();
}

PolyReducer::~PolyReducer()
{
  delete f;
  f = 0;
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////

void PolyReducer::insertTail(const_term multiplier, const Poly *g1)
{
  size_t ncmps = 0;

  if (g1->nTerms() <= 1) return;
  Poly *g = g1->copy();
  g->multByTerm(multiplier.coeff, multiplier.monom);

  Poly::iterator ig = g->begin();
  ++ig;
  Poly *h = Poly::add(R, f_iter, f->end(), ig, g->end(), ncmps);
  stats_n_compares += ncmps;
  stats_n_inserts++;

  delete f;
  delete g;
  f = h;
  f_iter = f->begin();

  size_t fmem = f->getMemoryUse();
  if (fmem > mMemUsage)
    mMemUsage = fmem;
}

void PolyReducer::insert(monomial multiplier, const Poly *g1)
{
  size_t ncmps = 0;

  Poly *g = g1->copy();
  g->multByMonomial(multiplier);

  Poly::iterator ig = g->begin();
  Poly *h = Poly::add(R, f_iter, f->end(), ig, g->end(), ncmps);
  stats_n_compares += ncmps;
  stats_n_inserts++;

  delete f;
  delete g;
  f = h;
  f_iter = f->begin();

  size_t fmem = f->getMemoryUse();
  if (fmem > mMemUsage)
    mMemUsage = fmem;
}

bool PolyReducer::leadTerm(const_term &result)
{
  if (f_iter != f->end())
    {
      result.coeff = f_iter.getCoefficient();
      result.monom = f_iter.getMonomial();
      return true;
    }
  return false;
}

void PolyReducer::removeLeadTerm()
// returns true if there is a term to extract
{
  if (f_iter == f->end()) return;
  ++f_iter;
}

void PolyReducer::value(Poly &result)
// keep extracting lead term until done
{
  Poly::iterator fend = f->end();
  result.append(f_iter, fend);
  delete f;
  f = new Poly(*R);
  f_iter = f->begin();
}

void PolyReducer::resetReducer()
{
  delete f;
  f = new Poly(*R);
  f_iter = f->begin();
}

void PolyReducer::dump() const
{
}

size_t PolyReducer::getMemoryUse() const {
  return mMemUsage;
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
