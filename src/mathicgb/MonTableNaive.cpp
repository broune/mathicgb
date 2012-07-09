// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include <ostream>
#include "MonTableNaive.hpp"
#include <iostream>

MonTableNaive::MonTableNaive(const PolyRing *R) :
  conf_(R),
  mPool(sizeof(mon_node)),
  table(0)
{
  stats_.n_member = 0;
  stats_.n_inserts = 0;
  stats_.n_insert_already_there = 0;
  stats_.n_compares = 0;
}

MonTableNaive::MonTableNaive(const Configuration& C) :
  conf_(C),
  mPool(sizeof(mon_node)),
  table(0)
{
  stats_.n_member = 0;
  stats_.n_inserts = 0;
  stats_.n_insert_already_there = 0;
  stats_.n_compares = 0;
}

MonTableNaive::~MonTableNaive()
{
  // Nothing to do?
}

bool MonTableNaive::member(const_monomial t, ValueType &result_val) const
{
  result_val = 0;
  stats_.n_member++;
  for (mon_node *p = table; p != 0; p = p->next)
    {
      stats_.n_compares++;
      switch (getPolyRing()->monomialCompare(p->monom, t))
        {
        case GT:
          return false;
        case EQ:
          // Already exists, leave
          return true;
        case LT:
          // Check divisibility
          if (getPolyRing()->monomialIsDivisibleBy(t, p->monom))
            return true;
        }
    }
  return false;
}

void MonTableNaive::insert_node(mon_node *p, const_monomial t)
{
  mon_node *q = static_cast<mon_node *>(mPool.alloc());
  q->next = p->next;
  p->next = q;
  q->monom = t;

  // Now remove the elements of the list p->next which are
  // divisible by q->monom

  mon_node *r = q;
  while (r->next != 0)
    {
      if (getPolyRing()->monomialIsDivisibleBy(r->next->monom, q->monom))
        {
          // remove this term:
          mon_node *a = r->next;
          r->next = a->next;
          mPool.free(a);
        }
      else
        r = r->next;
    }
}

bool MonTableNaive::insert(const_monomial t, ValueType /* val_ignored */)
{
  stats_.n_inserts++;
  // Only inserts if not in the table...
  mon_node head;
  head.next = table;
  mon_node *f = &head;
  while (f->next != 0)
    {
      stats_.n_compares++;
      switch (getPolyRing()->monomialCompare(f->next->monom, t))
        {
        case GT:
          // Insert t before f->next
          insert_node(f, t);
          table = head.next;
          return true;
        case EQ:
          // Already exists, leave
          stats_.n_insert_already_there++;
          return false;
        case LT:
          // Check divisibility
          if (getPolyRing()->monomialIsDivisibleBy(t, f->next->monom))
            {
              stats_.n_insert_already_there++;
              return false;
            }
          f = f->next;
        }
    }
  insert_node(f, t);
  table = head.next;
  return true;
}

size_t MonTableNaive::n_elems() const
{
  size_t len = 0;
  for (mon_node *q = table; q != 0; q = q->next) len++;
  return len;
}

void MonTableNaive::getMonomials(std::vector<const_monomial>& monomials) {
  for (mon_node *p = table; p != 0; p=p->next)
    monomials.push_back(p->monom);
}

void MonTableNaive::display(std::ostream &o, int /* level */) const
{
  o << n_elems() << ": ";
  for (mon_node *p = table; p != 0; p=p->next)
    {
      getPolyRing()->monomialDisplay(o, p->monom, false);
      o << "  ";
    }
  o << std::endl;
}

void MonTableNaive::displayStats(std::ostream &o) const
{
  o << "  #elements: " << n_elems() <<  std::endl;
  o << "  #calls member: " << stats_.n_member << std::endl;
  o << "  #calls insert, but already there: " << stats_.n_insert_already_there << std::endl;
  o << "  #compares: " << stats_.n_compares << std::endl;
}

void MonTableNaive::dump(int level) const
{
  displayStats(std::cerr);
  if (level > 0) display(std::cerr, level-1);
}

std::string MonTableNaive::getName() const
{
  return "naive";
}

size_t MonTableNaive::getMemoryUse() const
{
  return mPool.getMemoryUse();
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
