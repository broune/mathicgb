#ifndef _poly_basis_h_
#define _poly_basis_h_

#include "Poly.hpp"
#include "DivisorLookup.hpp"
#include <vector>
#include <memory>

class PolyRing;
class Ideal;
class FreeModuleOrder;

class PolyBasis {
public:
  // Ring and order must live for as long as this object. divisorLookupFactory
  // only needs to live for the duration of the constructor.
  PolyBasis(
    const PolyRing& ring,
    FreeModuleOrder& order,
    std::unique_ptr<DivisorLookup> divisorLookup);

  // Deletes the Poly's stored in the basis.
  ~PolyBasis();

  // Returns the initial monomial ideal of the basis (not the ideal).
  std::unique_ptr<Ideal> initialIdeal() const;

  // Inserts a polynomial into the basis at index size().
  // Lead monomials must be unique among basis elements.
  // So the index is size() - 1 afterwards since size() will increase by 1.
  void insert(std::unique_ptr<Poly> poly);

  // Returns the index of a basis element whose lead term divides mon.
  // Returns -1 if there is no such basis element.
  size_t divisor(const_monomial mon) const;

  // As the non-slow version, but uses simpler and slower code.
  size_t divisorSlow(const_monomial mon) const;

  // Replaces basis element at index with the given new value. The lead
  // term of the new polynomial must be the same as the previous one.
  // This is useful for auto-tail-reduction.
  void replaceSameLeadTerm(size_t index, std::unique_ptr<Poly> newValue) {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    MATHICGB_ASSERT(newValue.get() != 0);
    MATHICGB_ASSERT(!newValue->isZero());
    MATHICGB_ASSERT(mRing.monomialEQ
                    (leadMonomial(index), newValue->getLeadMonomial()));
    mDivisorLookup->remove(leadMonomial(index));
    delete mEntries[index].poly;
    mEntries[index].poly = newValue.release();
    mDivisorLookup->insert(leadMonomial(index), index);    
    MATHICGB_ASSERT(mEntries[index].poly != 0);
  }

  struct Stats {
    Stats():
      buchbergerLcmQueries(0),
      buchbergerLcmHits(0),
      buchbergerLcmNearHits(0),
      buchbergerLcmCacheHits(0) {}

    unsigned long long buchbergerLcmQueries;
    unsigned long long buchbergerLcmHits;
    unsigned long long buchbergerLcmNearHits;
    unsigned long long buchbergerLcmCacheHits;
  };
  Stats stats() const {return mStats;}

  // Returns true if Buchberger's second criterion for eliminating useless
  // S-pairs applies to the pair (a,b). Let
  //   l(a,b) = lcm(lead(a), lead(b)).
  // The criterion says that if there is some other basis element c such that
  //   lead(c)|l(a,b)
  // and
  //   l(a,c) reduces to zero, and
  //   l(b,c) reduces to zero
  // then (a,b) will reduce to zero (using classic non-signature reduction).
  //
  // This criterion is less straight forward to apply in case for example
  //   l(a,b) = l(a,c) = l(b,c)
  // since then there is the potential to erroneously eliminate all the three
  // pairs among a,b,c on the assumption that the other two pairs will reduce
  // to zero. In such cases, we eliminate the pair with the lowest indexes.
  // This allows removing generators that get non-minimal lead term without
  // problems.
  bool buchbergerLcmCriterion(size_t a, size_t b) const;

  // As the overload with an lcmAb parameter, except lcmAB must be the lcm of
  // the lead monomials of a and b. Then this quantity does not have to be
  // computed for the criterion.
  bool buchbergerLcmCriterion(size_t a, size_t b, const_monomial lcmAB) const;

  // As the non-slow version, but uses simpler and slower code.
  bool buchbergerLcmCriterionSlow(size_t a, size_t b) const;

  // Returns the number of basis elements, including retired elements.
  size_t size() const {return mEntries.size();}

  // Returns the ambient polynomial ring of the polynomials in the basis.
  const PolyRing& ring() const {return mRing;}

  // Returns the term order on the basis.
  const FreeModuleOrder& order() const {return mOrder;}

  // Returns a data structure containing the lead monomial of each lead
  // monomial.
  const DivisorLookup& divisorLookup() const {return *mDivisorLookup;}

  // Retires the basis element at index, which frees the memory associated
  // to it, including the basis element polynomial, and marks it as retired. 
  // todo: implement
  std::unique_ptr<Poly> retire(size_t index) {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    mDivisorLookup->remove(leadMonomial(index));
    std::unique_ptr<Poly> poly(mEntries[index].poly);
    mEntries[index].poly = 0;
    mEntries[index].retired = true;
    return poly;
  }

  // Returns true of the basis element at index has been retired.
  bool retired(size_t index) const {
    MATHICGB_ASSERT(index < size());
    return mEntries[index].retired;
  }

  // Returns the basis element polynomial at index.
  Poly& poly(size_t index) {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    return *mEntries[index].poly;
  }

  // Returns the basis element polynomial at index.
  const Poly& poly(size_t index) const {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    return *mEntries[index].poly;
  }

  const_monomial leadMonomial(size_t index) const {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    return poly(index).getLeadMonomial();
  }

  coefficient leadCoefficient(size_t index) const {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    return poly(index).getLeadCoefficient();
  }

  // Returns true if the leading monomial of the basis element at index is not
  // divisible by the lead monomial of any other basis element. Lead
  // monomials are required to be unique among basis elements, so the case
  // of several equal lead monomials does not occur.
  bool leadMinimal(size_t index) const {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    MATHICGB_SLOW_ASSERT(mEntries[index].leadMinimal == leadMinimalSlow(index));
    return mEntries[index].leadMinimal;
  }

  // Returns true if m is not divisible by the lead monomial of any
  // basis element. Equality counts as divisibility.
  bool leadMinimal(const Poly& poly) const {
    MATHICGB_ASSERT(&poly != 0);
    return mDivisorLookup->divisor(poly.getLeadMonomial()) !=
      static_cast<size_t>(-1);
  }

  // Returns the number of basis elements with minimal lead monomial.
  size_t minimalLeadCount() const;

  // Returns the index of the basis element of maximal index
  // whose lead monomial is minimal.
  size_t maxIndexMinimalLead() const;

  // Returns the basis element polynomial at index.
  const Poly& basisElement(size_t index) const {
    MATHICGB_ASSERT(index < size());
    MATHICGB_ASSERT(!retired(index));
    return *mEntries[index].poly;
  }

  // Returns the number of monomials across all the basis elements.
  // Monomials that appear in more than one basis element are counted more
  // than once.
  size_t monomialCount() const;

  // Returns how many bytes has been allocated by this object.
  size_t getMemoryUse() const;

  void usedAsStart(size_t index) const {
    MATHICGB_ASSERT(index < size());
    ++mEntries[index].usedAsStartCount;
  }

  unsigned long long usedAsStartCount(size_t index) const {
    MATHICGB_ASSERT(index < size());
    return mEntries[index].usedAsStartCount;
  }

  void usedAsReducer(size_t index) const {
    MATHICGB_ASSERT(index < size());
    ++mEntries[index].usedAsReducerCount;
  }

  unsigned long long usedAsReducerCount(size_t index) const {
    MATHICGB_ASSERT(index < size());
    return mEntries[index].usedAsReducerCount;
  }

  void wasPossibleReducer(size_t index) const {
    MATHICGB_ASSERT(index < size());
    ++mEntries[index].possibleReducerCount;
  }

  unsigned long long wasPossibleReducerCount(size_t index) const {
    MATHICGB_ASSERT(index < size());
    return mEntries[index].possibleReducerCount;
  }

  void wasNonSignatureReducer(size_t index) const {
    MATHICGB_ASSERT(index < size());
    ++mEntries[index].nonSignatureReducerCount;
  }

  unsigned long long wasNonSignatureReducerCount(size_t index) const {
    MATHICGB_ASSERT(index < size());
    return mEntries[index].nonSignatureReducerCount;
  }

private:
  // Slow versions use simpler code. Used to check results in debug mode.
  bool leadMinimalSlow(size_t index) const;

  class Entry {
  public:
    Entry();

    Poly* poly;
    bool leadMinimal;
    bool retired;

    // Statistics on reducer choice in reduction
    mutable unsigned long long usedAsStartCount;
    mutable unsigned long long usedAsReducerCount;
    mutable unsigned long long possibleReducerCount;
    mutable unsigned long long nonSignatureReducerCount;
  };
  typedef std::vector<Entry> EntryCont;
  typedef EntryCont::iterator EntryIter;
  typedef EntryCont::const_iterator EntryCIter;

  const PolyRing& mRing;
  FreeModuleOrder& mOrder;
  std::unique_ptr<DivisorLookup> mDivisorLookup;
  std::vector<Entry> mEntries;
  mutable Stats mStats;

  static const bool mUseBuchbergerLcmHitCache = true;
  mutable std::vector<size_t> mBuchbergerLcmHitCache;
};

#endif
