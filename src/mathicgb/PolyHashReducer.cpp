// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "PolyHashReducer.hpp"

MATHICGB_NAMESPACE_BEGIN

PolyHashReducer::PolyHashReducer(const PolyRing *R0):
  R_(R0), H_(R0,15)
{
  f_ = new HashPoly;
  f_iter_ = f_->begin();
}

PolyHashReducer::~PolyHashReducer()
{
  delete f_;
  f_ = 0;
}

void PolyHashReducer::merge(const HashPoly::const_iterator &fbegin, const HashPoly::const_iterator &fend,
                            const HashPoly::const_iterator &gbegin, const HashPoly::const_iterator &gend,
                            HashPoly &result,
                            size_t &n_compares)
{
  HashPoly::const_iterator f = fbegin;
  HashPoly::const_iterator g = gbegin;

  n_compares = 0;

  if (f == fend)
    result.insert(result.end(), g, gend);
  else if (g == gend)
    result.insert(result.end(), f, fend);
  else {
    bool done = false;
    while (!done) {
      int cmp = R_->monomialCompare((*f)->monom, (*g)->monom);
      n_compares++;
      MATHICGB_ASSERT(cmp != EQ);
      if (cmp == LT) {
        result.push_back(*g);
        ++g;
        if (g == gend)
          {
            result.insert(result.end(), f, fend);
            done = true;
          }
      } else {
        MATHICGB_ASSERT(cmp == GT);
        result.push_back(*f);
        ++f;
        if (f == fend) {
          result.insert(result.end(), g, gend);
          done = true;
        }
      }
    }
  }
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////

void PolyHashReducer::insertTail(const_term multiplier, const Poly *g1)
{
  size_t ncmps = 0;

  if (g1->nTerms() <= 1) return;

  HashPoly M;
  //  Poly *g = g1->copy(); //MES: ONLY COPY THE TAIL!!
  //  g->multByTerm(multiplier.coeff, multiplier.monom);
  //  H_.fromPoly(*g, M);

  H_.insert(multiplier, ++(g1->begin()), g1->end(), M);

  if (!M.empty())
    {
      HashPoly h;
      merge(f_iter_, f_->end(), M.begin(), M.end(), h, ncmps);
      swap(*f_, h);
      f_iter_ = f_->begin();
    }

  stats_n_compares += ncmps;
  stats_n_inserts++;
  // delete g;
}

void PolyHashReducer::insert(monomial multiplier, const Poly *g1)
{
  size_t ncmps = 0;

  //  Poly *g = g1->copy();
  //  g->multByMonomial(multiplier);

  HashPoly M;
  //  H_.fromPoly(*g, M);

  H_.insert(multiplier, g1->begin(), g1->end(), M);

  HashPoly h;
  merge(f_iter_, f_->end(), M.begin(), M.end(), h, ncmps);
  stats_n_compares += ncmps;
  stats_n_inserts++;

  //  delete g;
  swap(*f_, h);
  f_iter_ = f_->begin();
}

bool PolyHashReducer::leadTerm(const_term &result)
{
  while (f_iter_ != f_->end())
    {
      if (H_.popTerm(*f_iter_, result.coeff, result.monom))
        // returns true *f_iter is not the zero element
        return true;
      ++f_iter_;
    }
  return false;
}

void PolyHashReducer::removeLeadTerm()
// returns true if there is a term to extract
{
  if (f_iter_ == f_->end()) return;
  ++f_iter_;
}

void PolyHashReducer::value(Poly &result)
// keep extracting lead term until done
{
  const_term t;
  for ( ; f_iter_ != f_->end(); ++f_iter_)
    if (leadTerm(t))
      result.appendTerm(t.coeff, t.monom);
  resetReducer();
}

size_t PolyHashReducer::getMemoryUse() const
{
  return
    TypicalReducer::getMemoryUse() +
    H_.getMemoryUse() +
    f_->capacity() * sizeof(PolyHashTable::node *);
}

void PolyHashReducer::resetReducer()
{
  const_term t;
  for ( ; f_iter_ != f_->end(); ++f_iter_)
    leadTerm(t);

  delete f_;
  f_ = new HashPoly;
  f_iter_ = f_->begin();

  H_.reset();
}

void PolyHashReducer::dump() const
{
  H_.dump(0);
}

MATHICGB_NAMESPACE_END
