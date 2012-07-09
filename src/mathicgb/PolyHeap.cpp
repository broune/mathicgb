// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include "PolyHeap.hpp"

#include <iostream>
#include <algorithm>

extern int tracingLevel;

size_t PolyHeap::stats_static_n_compares = 0;

PolyHeap::PolyHeap(const PolyRing *R_)
  : Reducer(),
    R(R_),
    heapCompare(R_),
    lead_is_computed(false)
{
  stats_static_n_compares = 0;
}

void PolyHeap::insert(const_term multiplier, Poly::const_iterator first, Poly::const_iterator last)
{
  if (first == last) return;
  heap_term t;
  t.multiplier = multiplier;
  t.first = first;
  t.last = last;
  t.actual = R->allocMonomial(mArena);
  R->monomialMult(multiplier.monom, t.first.getMonomial(), t.actual);
  terms_.push_back(t);
  std::push_heap(terms_.begin(), terms_.end(), heapCompare);
  stats_n_compares += stats_static_n_compares;
  stats_static_n_compares = 0;
  if (terms_.size() > stats_maxsize)
    stats_maxsize = terms_.size();
  stats_n_inserts++;
}

void PolyHeap::insertTail(const_term multiplier, const Poly *f)
{
  if (tracingLevel > 100) {
    std::cerr << "PolyHeap inserting tail of ";
    f->display(std::cerr);
    std::cerr << std::endl;
    std::cerr << "multiplied by: " << multiplier.coeff << "  *  ";
    R->monomialDisplay(std::cerr, multiplier.monom);
    std::cerr << std::endl;
  }

  Poly::const_iterator a = f->begin();
  insert(multiplier, ++a, f->end());
}
void PolyHeap::insert(monomial multiplier, const Poly *f)
{
  const_term t;
  t.monom = multiplier;
  R->coefficientSetOne(t.coeff);
  Poly::const_iterator a = f->begin();
  insert(t, a, f->end());
}

void PolyHeap::extract1(const_term &t)
{
  // Pop the largest term into t, incrememnt that poly iterator
  std::pop_heap(terms_.begin(), terms_.end(), heapCompare);
  stats_n_compares += stats_static_n_compares;
  stats_static_n_compares = 0;
  heap_term h = *(terms_.rbegin());
  t.monom = h.actual;
  t.coeff = h.multiplier.coeff;
  R->coefficientMultTo(t.coeff, h.first.getCoefficient());
  terms_.pop_back(); // removes h from heap
  ++h.first;
  if (h.first != h.last)
    insert(h.multiplier, h.first, h.last);
  else {
    // REMOVE h.multiplier
  }
}

bool PolyHeap::extractLeadTerm(const_term &result)
{
  bool result_set = false;
  for (;;) {
    if (terms_.empty()) break;
    if (!result_set)
      {
        extract1(result);
        result_set = true;
        continue;
      }
    heap_term &h = terms_[0]; // Now look at the next term
    if (!R->monomialEQ(h.actual, result.monom)) break;
    // at this point, we need to grab the monomial
    const_term t;
    extract1(t);
    R->coefficientAddTo(result.coeff, t.coeff);
    result_set = !(R->coefficientIsZero(result.coeff));
  }
  return result_set;
}

bool PolyHeap::computeLeadTerm()
{
  if (!lead_is_computed)
    {
      for (;;) {
        if (terms_.empty()) break;
        if (!lead_is_computed)
          {
            extract1(lead_term);
            lead_is_computed = true;
            continue;
          }
        heap_term &h = terms_[0]; // Now look at the next term
        if (!R->monomialEQ(h.actual, lead_term.monom)) break;
        // at this point, we need to grab the monomial
        const_term t;
        extract1(t);
        R->coefficientAddTo(lead_term.coeff, t.coeff);
        lead_is_computed = !(R->coefficientIsZero(lead_term.coeff));
      }
    }
  return lead_is_computed; // will be false if could not compute
}

bool PolyHeap::findLeadTerm(const_term &result)
{
  if (!lead_is_computed && !computeLeadTerm()) return false;
  result = lead_term;
  return true;
}

void PolyHeap::removeLeadTerm()
{
  if (!lead_is_computed) computeLeadTerm();
  lead_is_computed = false;
}

void PolyHeap::value(Poly &result)
{
  const_term t;
  while (findLeadTerm(t))
    {
      result.appendTerm(t.coeff, t.monom);
      lead_is_computed = false;
    }
  resetReducer();
}

size_t PolyHeap::getMemoryUse() const
{
  return
    Reducer::getMemoryUse() +
    terms_.capacity() * sizeof(heap_term);
}

void PolyHeap::resetReducer()
{
  terms_.resize(0);
  lead_is_computed = false;
}

void PolyHeap::dump() const
// display debugging info about the heap
{
  std::cout << "-- polyheap --" << std::endl;
  for (size_t i = 0; i < terms_.size(); i++)
    {
      std::cout << "actual: ";
      R->monomialDisplay(std::cout, terms_[i].actual);

      std::cout << " monom: ";
      R->monomialDisplay(std::cout, terms_[i].multiplier.monom);

      std::cout << " coeff: " << terms_[i].multiplier.coeff ;
      std::cout << std::endl;
    }
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
