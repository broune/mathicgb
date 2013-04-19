// Copyright 2011 Michael E. Stillman

#ifndef _groebner_basis_h_
#define _groebner_basis_h_

#include <vector>
#include <set>
#include "PolyRing.hpp"
#include "Poly.hpp"
#include "FreeModuleOrder.hpp"
#include "DivisorLookup.hpp"
#include "PolyBasis.hpp"
#include "MonoProcessor.hpp"

#ifndef USE_RATIO_RANK
#define USE_RATIO_RANK true
#endif

class GroebnerBasis {
public:
  GroebnerBasis(
    const PolyRing* R,
    FreeModuleOrder* order,
    int divisorLookupType,
    int monTableType,
    bool preferSparseReducers
  );
  ~GroebnerBasis();

  const PolyRing& ring() const {return mBasis.ring();}
  const Poly& poly(size_t index) const {return mBasis.poly(index);}
  size_t size() const {return mBasis.size();}
  const FreeModuleOrder& order() const {return mBasis.order();}

  // todo: stop using non-const version of basis() and then remove it.
  PolyBasis& basis() {return mBasis;}
  const PolyBasis& basis() const {return mBasis;}

  const_monomial getLeadMonomial(size_t i) const {
    return mBasis.leadMonomial(i);
  }

  coefficient getLeadCoefficient(size_t i) const {
    return mBasis.leadCoefficient(i);
  }

  const_monomial getSigLeadRatio(size_t i) const {
    MATHICGB_ASSERT(i < size());
    return sigLeadRatio[i];
  }

  // Signifies that the module has taken on another e_i.
  // Must call this before adding a polynomial to the basis with
  // a signature in the new component.
  void addComponent();

  // Takes over ownership of sig and f. sig must come from the pool
  // of ring() and f must have been allocated with new.
  void insert(monomial sig, std::unique_ptr<Poly> f);

  const_monomial getSignature(size_t i) const {
    MATHICGB_ASSERT(i < size());
    return mSignatures[i];
  }

  // Returns the index of a basis element that regular reduces term in
  // signature sig. Returns -1 if no such element exists. A basis element
  // u is a regular reducer if leadTerm(u) divides term
  // and (term / leadTerm(u)) * signature(u) < sig.
  size_t regularReducer(const_monomial sig, const_monomial term) const;

  // Uses the functionality in the divisor finder for
  // computing up to maxDivisors low ratio base divisors.
  // The divisors are placed into divisors.
  void lowBaseDivisors(
    std::vector<size_t>& divisors,
    size_t maxDivisors,
    size_t newGenerator) const;

  // Uses the functionality in the divisor finder for
  // computing a high base divisor. Returns the index
  // of the divisor or -1 if none are found.
  size_t highBaseDivisor(size_t newGenerator) const;

  // Find the basis element g_i whose signature S divides sig
  // such that (S/sig)g_i has minimal leading term. Returns i.
  size_t minimalLeadInSig(const_monomial sig) const;

  // Returns true if poly can be singular reduced in signature sig.
  // In other words, returns true if there is a basis element with
  // lead term M and signature S such that M divides the lead term N
  // of poly and such that N/M*S == sig. This is slow because it is
  // currently only used for asserts - implement a fast version if
  // that changes.
  bool isSingularTopReducibleSlow(const Poly& poly, const_monomial sig) const;

  void display(std::ostream &o) const;
  void displayBrief(std::ostream &o) const;
  void dump() const;
  size_t getMemoryUse() const;

  // Compares the signature/lead ratio of basis element a to basis element b
  // and returns LT, EQ or GT.
  inline int ratioCompare(size_t a, size_t b) const;

  /// Post processes all signatures. This currently breaks all sorts
  /// of internal invariants - it's supposed to be a temporary hack.
  void postprocess(const MonoProcessor<PolyRing::Monoid>& processor);

  class StoredRatioCmp {
  public:
    // Stores the ratio numerator/denominator and prepares it for comparing
    // to the sig/lead ratios in basis.
    StoredRatioCmp(
      const_monomial numerator,
      const_monomial denominator,
      const GroebnerBasis& basis);
    ~StoredRatioCmp();

    // compares the stored ratio to the basis element with index be.
    inline int compare(size_t be) const;

  private:
    StoredRatioCmp(const StoredRatioCmp&); // not available
    void operator=(const StoredRatioCmp&); // not available

    const GroebnerBasis& mBasis;
    size_t mRatioRank;
    monomial mRatio;
    mutable monomial mTmp;
  };

private:
  // Slow versions use simpler code. Used to check results in debug mode.
  size_t regularReducerSlow(const_monomial sig, const_monomial term) const;
  size_t minimalLeadInSigSlow(const_monomial sig) const;
  size_t highBaseDivisorSlow(size_t newGenerator) const;
  void lowBaseDivisorsSlow(
    std::vector<size_t>& divisors,
    size_t maxDivisors,
    size_t newGenerator) const;

  friend class StoredRatioCmp;

  const DivisorLookup& divisorLookup() const {
    return mBasis.divisorLookup();
  }

  std::unique_ptr<DivisorLookup::Factory const> const mDivisorLookupFactory;

  // may change at next insert!
  size_t ratioRank(size_t index) const {
    MATHICGB_ASSERT(index < size());
    return mRatioRanks[index];
  }

  // Only useful for comparing to basis elements. Two ratios might get the same
  // rank without being equal. The proper rank may change when a new generator
  // is added.
  size_t ratioRank(const_monomial ratio) const;

  std::vector<monomial> mSignatures;

  // the ratio signature/initial term including negative entries and module component
  std::vector<monomial> sigLeadRatio;

  // true if giving each generator an integer id based on its
  // position in a sorted order of sig/lead ratios.
  static const bool mUseRatioRank = USE_RATIO_RANK;
  static const bool mUseStoredRatioRank = USE_RATIO_RANK;
  class RatioOrder {
  public:
    RatioOrder(std::vector<monomial>& ratio, const FreeModuleOrder& order):
        mRatio(ratio), mOrder(order) {}
    bool operator()(size_t a, size_t b) const {
      return mOrder.signatureCompare(mRatio[a], mRatio[b]) == LT;
    }
  private:
    std::vector<monomial>& mRatio;
    const FreeModuleOrder& mOrder;
  };
  typedef std::multiset<size_t, RatioOrder> RatioSortedType;
  typedef size_t Rank;
  RatioSortedType mRatioSorted;
  std::vector<Rank> mRatioRanks;

  std::vector<DivisorLookup*> mSignatureLookup;

  // Contains those lead terms that are minimal.
  std::unique_ptr<DivisorLookup> const mMinimalDivisorLookup;

  PolyBasis mBasis;
  bool const mPreferSparseReducers;
  mutable monomial mTmp;
};

inline int GroebnerBasis::ratioCompare(size_t a, size_t b) const {
  if (mUseRatioRank) {
#ifdef MATHICGB_DEBUG
    int const value =
      order().signatureCompare(getSigLeadRatio(a), getSigLeadRatio(b));
#endif
    if (mRatioRanks[a] < mRatioRanks[b]) {
      MATHICGB_ASSERT_NO_ASSUME(value == LT);
      return LT;
    } else if (mRatioRanks[a] > mRatioRanks[b]) {
      MATHICGB_ASSERT_NO_ASSUME(value == GT);
      return GT;
    } else {
      MATHICGB_ASSERT_NO_ASSUME(value == EQ);
      return EQ;
    }
  } else {
    // A/a < B/b   <=>  A < (B/b)a
    ring().monomialDivideToNegative(getSignature(b), getLeadMonomial(b), mTmp);
    ring().monomialMultTo(mTmp, getLeadMonomial(a));
    int value = order().signatureCompare(getSignature(a), mTmp);
    MATHICGB_ASSERT(value ==
      order().signatureCompare(getSigLeadRatio(a), getSigLeadRatio(b)));
    return value;
  }
}

inline int GroebnerBasis::StoredRatioCmp::compare(size_t be) const {
  if (GroebnerBasis::mUseStoredRatioRank) {
#ifdef DEBUG
    int value = mBasis.order().
      signatureCompare(mRatio, mBasis.getSigLeadRatio(be));
#endif
    GroebnerBasis::Rank otherRank = mBasis.ratioRank(be);
    if (mRatioRank < otherRank) {
      MATHICGB_ASSERT_NO_ASSUME(value == LT);
      return LT;
    } else if (mRatioRank > otherRank) {
      MATHICGB_ASSERT_NO_ASSUME(value == GT);
      return GT;
    } else {
      MATHICGB_ASSERT_NO_ASSUME(value == EQ);
      return EQ;
    }
  } else {
    mBasis.ring().monomialMult(mRatio, mBasis.getLeadMonomial(be), mTmp);
    int value = mBasis.order().signatureCompare(mTmp, mBasis.getSignature(be));
    MATHICGB_ASSERT(value ==
      mBasis.order().signatureCompare(mRatio, mBasis.getSigLeadRatio(be)));
    return value;
  }
}

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
