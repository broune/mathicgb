// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_TOURNAMENT_REDUCER_GUARD
#define MATHICGB_TOURNAMENT_REDUCER_GUARD

#include "TypicalReducer.hpp"
#include <mathic.h>
#include <memtailor.h>

MATHICGB_NAMESPACE_BEGIN

class TournamentReducer : public TypicalReducer {
public:
  TournamentReducer(const PolyRing& R);
  virtual ~TournamentReducer();

  virtual std::string description() const { return "tournament tree reducer"; }

  virtual void insertTail(const_term multiplier, const Poly *f);
  virtual void insert(monomial multiplier, const Poly *f);

  virtual bool leadTerm(const_term &result);
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

    void computeCurrent(const PolyRing& ring);
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
  
  void insert(const_term multiplier, Poly::iterator first, Poly::iterator last);

  const PolyRing& mRing;
  term mLeadTerm;
  bool mLeadTermKnown;
  mic::TourTree<Configuration> mQueue;
  memt::BufferPool mPool;
};

MATHICGB_NAMESPACE_END
#endif
