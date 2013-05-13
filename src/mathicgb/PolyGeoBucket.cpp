// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "stdinc.h"
#include "PolyGeoBucket.hpp"

#include <mathic.h>

MATHICGB_NAMESPACE_BEGIN

const size_t heap_size[GEOHEAP_SIZE] = {4, 16, 64, 256, 1024, 4096,
                               16384, 65536, 262144, 1048576, 4194304,
                               16777216, 67108864, 268435456,
                               1073741824};

PolyGeoBucket::PolyGeoBucket(const PolyRing *R0): 
  R(R0), top_of_heap(-1), lead(-1)
{
  int i;
  for (i=0; i<GEOHEAP_SIZE; i++)
    {
      heap[i].poly = 0;
      // the iterators are NOT VALID if this poly is the null value
    }
}

PolyGeoBucket::~PolyGeoBucket()
{
  resetReducer();
}

///////////////////////////////////////
// Helper routines ////////////////////
///////////////////////////////////////

void PolyGeoBucket::add_to(heap_record &a, heap_record &b)
// handles the case when a.poly == 0 too
// If the result is 0, then a.poly is set to 0.
{
  assert(b.poly != 0);
  assert(b.first != b.last);
  if (a.poly == 0)
    {
      a = b;
      b.poly = 0;
    }
  else
    {
      size_t ncmps;
      Poly *g = Poly::add(R, a.first, a.last, b.first, b.last, ncmps);
      stats_n_inserts++;
      delete a.poly;
      a.poly = 0;
      stats_n_compares += ncmps;
      if (g->isZero())
        {
          delete g;
        }
      else
        {
          a.poly = g;
          a.first = g->begin();
          a.last = g->end();
        }
      delete b.poly;
      b.poly = 0;
    }
}

void PolyGeoBucket::update_size_stats()
{
  size_t size_polys = 0;
  size_t size_live = 0;
  for (int i = 0; i<= top_of_heap; i++)
    {
      if (heap[i].poly == 0) continue;
      size_polys += heap[i].poly->nTerms();
      size_live += heap[i].last - heap[i].first;
    }
  if (size_polys > stats_maxsize)
    stats_maxsize = size_polys;
  if(size_live >= stats_maxsize_live)
    stats_maxsize_live = size_live;
}
void PolyGeoBucket::insert(heap_record &a)
{
  // Takes the heap_record, and adds it into one of the heap elems
  // If the size changes, then this continues until the heap
  //  sizes are valid

  lead = -1;
  size_t len = a.last - a.first;
  int i= 0;
  while (len >= heap_size[i]) i++;

  add_to(heap[i], a);
  if (heap[i].poly != 0)
    {
      len = heap[i].last - heap[i].first;
      while (len >= heap_size[i])
        {
          i++;

          add_to(heap[i], heap[i-1]);
          if (heap[i].poly == 0) break;
          len = heap[i].last - heap[i].first;
        }
    }
  if (i > top_of_heap) {
    top_of_heap = i;
    if (i >= GEOHEAP_SIZE)
      mic::reportInternalError("Too many levels in PolyGeoBucket.");
  }

  stats_n_inserts++;
  update_size_stats();
}

bool PolyGeoBucket::computeLeadTerm()
// returns true if the heap does not add to 0
{
  int lead_so_far = -1;
  for (int i=0; i <= top_of_heap; i++)
    {
      if (heap[i].poly == 0) continue;
      if (lead_so_far < 0)
        {
          lead_so_far = i;
          continue;
        }
      int cmp = R->monomialCompare(heap[lead_so_far].first.getMonomial(), heap[i].first.getMonomial());
      stats_n_compares++;
      if (cmp == GT) continue;
      if (cmp == LT)
        {
          lead_so_far = i;
          continue;
        }
      // At this point we have equality
      R->coefficientAddTo(heap[lead_so_far].first.getCoefficient(), heap[i].first.getCoefficient());
      // now increment one of these
      ++heap[i].first;
      if (heap[i].first == heap[i].last)
        {
          delete heap[i].poly;
          heap[i].poly = 0;
        }
      if (R->coefficientIsZero(heap[lead_so_far].first.getCoefficient()))
        {
          ++heap[lead_so_far].first;
          if (heap[lead_so_far].first == heap[lead_so_far].last)
            {
              delete heap[lead_so_far].poly;
              heap[lead_so_far].poly = 0;
            }
          lead_so_far = -1;
          i = -1;
        }
    }
  lead = lead_so_far;
  return (lead_so_far >= 0);
}

///////////////////////////////////////
// External interface routines ////////
///////////////////////////////////////

void PolyGeoBucket::insertTail(const_term multiplier, const Poly *f)
{
  if (f->nTerms() <= 1) return;
  Poly *g = f->copy();
  g->multByTerm(multiplier.coeff, multiplier.monom);
  heap_record a;
  a.poly = g;
  a.first = g->begin();
  ++a.first;
  a.last = g->end();
  insert(a);
}

void PolyGeoBucket::insert(monomial multiplier, const Poly *f)
{
  Poly *g = f->copy();
  g->multByMonomial(multiplier);
  heap_record a;
  a.poly = g;
  a.first = g->begin();
  a.last = g->end();
  insert(a);
}

bool PolyGeoBucket::leadTerm(const_term &result)
{
  if (lead == -1 && !computeLeadTerm()) return false;
  result.coeff = heap[lead].first.getCoefficient();
  result.monom = heap[lead].first.getMonomial();
  return true;
}

void PolyGeoBucket::removeLeadTerm()
// returns true if there is a term to extract
{
  if (lead == -1 && !computeLeadTerm()) return;
  ++heap[lead].first;
  if (heap[lead].first == heap[lead].last)
    {
      // We are at the end here, so free poly, and set poly to 0.
      delete heap[lead].poly;
      heap[lead].poly = 0;
      lead = -1;
      return;
    }
  lead = -1;
}

void PolyGeoBucket::value(Poly &result)
// keep extracting lead term until done
{
  heap_record a;
  a.poly = 0;
  for (int i=0; i<=top_of_heap; i++)
    {
      if (heap[i].poly == 0) continue;
      add_to(a, heap[i]);
    }
  top_of_heap = -1;
  if (a.poly != 0)
    {
      result.append(a.first, a.last);
      delete a.poly;
      a.poly = 0;
    }
}

void PolyGeoBucket::resetReducer()
{
  for (int i=0; i<=top_of_heap; i++)
    {
      if (heap[i].poly != 0)
        {
          delete heap[i].poly;
          heap[i].poly = 0;
        }
    }
  top_of_heap = -1;
  lead = -1;
}

size_t PolyGeoBucket::getMemoryUse() const
{
  // size includes: each poly in the heap, as well as the
  // size of the heap itself
  size_t result =
    TypicalReducer::getMemoryUse() +
    sizeof(heap_record) * GEOHEAP_SIZE;
  for (int i=0; i<GEOHEAP_SIZE; ++i)
    if (heap[i].poly != 0)
      result += heap[i].poly->getMemoryUse();
  return result;
}

void PolyGeoBucket::dump() const
{
}

MATHICGB_NAMESPACE_END
