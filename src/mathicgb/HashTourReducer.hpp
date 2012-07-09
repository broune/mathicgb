// Copyright 2011 Bjarke Roune, Michael E. Stillman

#ifndef _hash_tour_reducer_h_
#define _hash_tour_reducer_h_

#include "Reducer.hpp"
#include "PolyHashTable.hpp"
#include <mathic.h>
#include <memtailor.h>

class HashTourReducer : public Reducer {
public:
  HashTourReducer(const PolyRing& R);
  virtual ~HashTourReducer();

  virtual std::string description() const {
    return "hashed tournament tree reducer";
  }

  virtual void insertTail(const_term multiplier, const Poly *f);
  virtual void insert(monomial multiplier, const Poly *f);

  virtual bool findLeadTerm(const_term &result);
  virtual void removeLeadTerm();

  virtual void value(Poly &result); // keep extracting lead term until done

  virtual size_t getMemoryUse() const;

  virtual void dump() const; // Used for debugging

protected:
  virtual void resetReducer();

private:
  // Represents a term multiple of a polynomial, 
  // together with a current term of the multiple.
  struct MultipleWithPos {
    MultipleWithPos(const Poly& poly, const_term multiple);

    Poly::const_iterator pos;
    Poly::const_iterator const end;
    const_term const multiple;
    monomial current; // multiple.monom * pos.getMonomial()
    PolyHashTable::node* node;

    void computeCurrent(const PolyRing& ring, monomial current);
    void currentCoefficient(const PolyRing& ring, coefficient& coeff);
    void destroy(const PolyRing& ring);
  };

  class Configuration {
  public:
    typedef MultipleWithPos* Entry;

    Configuration(const PolyRing& ring) : mRing(ring) {}

    typedef bool CompareResult;
    CompareResult compare(const Entry& a, const Entry& b) const {
      return mRing.monomialLT(a->current, b->current);
    }
    bool cmpLessThan(CompareResult r) const {return r;}

    static const bool fastIndex = true;

  private:
    const PolyRing& mRing;
  };

  class MonomialFree;

  void insertEntry(MultipleWithPos* entry);
  void insert(const_term multiplier, Poly::iterator first, Poly::iterator last);

  const PolyRing& mRing;
  term mLeadTerm;
  bool mLeadTermKnown;
  mic::TourTree<Configuration> mQueue;
  PolyHashTable mHashTable;
  memt::BufferPool mPool;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
