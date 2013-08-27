// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_DIV_LOOKUP_GUARD
#define MATHICGB_DIV_LOOKUP_GUARD

#include "SigPolyBasis.hpp"
#include "DivisorLookup.hpp"
#include "PolyRing.hpp"
#include <string>
#include <vector>
#include <iostream>

MATHICGB_NAMESPACE_BEGIN

/** Configuration class for interface to KDTree, DivList */
/* As such, it has entries that both will expect */
/* It also contains enough for the Naive monomial table use */

template<bool AllowRemovals, bool UseDivMask>
class DivLookupConfiguration;

template<bool AR, bool DM>
class DivLookupConfiguration {
public:
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::ConstMonoPtr ConstMonoPtr;
  typedef Monoid::ConstMonoRef ConstMonoRef;

  DivLookupConfiguration(
    const Monoid& monoid,
    int type,
    bool preferSparseReducers
  ):
    mBasis(0),
    mSigBasis(0),
    mMonoid(monoid),
    mType(type),
    mPreferSparseReducers(preferSparseReducers)
  {}

  void setBasis(const PolyBasis& basis) {
    if (mBasis == &basis)
      return;
    MATHICGB_ASSERT(mBasis == 0);
    MATHICGB_ASSERT(monoid() == basis.ring().monoid());
    mBasis = &basis;
  }

  void setSigBasis(const SigPolyBasis& sigBasis) {
    if (mSigBasis == &sigBasis)
      return;
    MATHICGB_ASSERT(mSigBasis == 0);
    MATHICGB_ASSERT(mBasis == 0 || mBasis == &sigBasis.basis());
    MATHICGB_ASSERT(monoid() == sigBasis.basis().ring().monoid());
    mSigBasis = &sigBasis;
    setBasis(sigBasis.basis());
  }

  const SigPolyBasis& sigBasis() const {
    MATHICGB_ASSERT(mSigBasis != 0);
    return *mSigBasis;
  }
  
  const PolyBasis& basis() const {
    MATHICGB_ASSERT(mBasis != 0);
    return *mBasis;
  }

  const Monoid& monoid() const {return mMonoid;}
  int type() const {return mType;}
  bool preferSparseReducers() const {return mPreferSparseReducers;}

  // *** Interface expected by mathic follows below

  typedef int Exponent;
  typedef ConstMonoRef Monomial;
  struct Entry {
    Entry(): monom(), index(static_cast<size_t>(-1)) {}
    Entry(ConstMonoRef monom0, size_t index0):
      monom(monom0.ptr()), index(index0) {}

    ConstMonoPtr monom;
    size_t index;
  };

  Exponent getExponent(const Monomial& m, size_t var) const {
    return monoid().exponent(m, var);
  }

  Exponent getExponent(const Entry& e, size_t var) const {
    return monoid().exponent(*e.monom, var);
  }

  bool divides(const Monomial& a, const Monomial& b) const {
    return monoid().divides(a, b);
  }

  bool divides(const Entry& a, const Monomial& b) const {
    return monoid().divides(*a.monom, b);
  }

  bool divides(const Monomial& a, const Entry& b) const {
    return monoid().divides(a, *b.monom);
  }

  bool divides(const Entry& a, const Entry& b) const {
    return monoid().divides(*a.monom, *b.monom);
  }

  bool getSortOnInsert() const {return false;}
  template<class A, class B>
  bool isLessThan(const A& a, const B& b) const {
    MATHICGB_ASSERT(false);
    return false;
  }

  size_t getVarCount() const {return monoid().varCount();}

  static const bool UseTreeDivMask = DM;
  static const bool UseLinkedList = false;
  static const bool UseDivMask = DM;
  static const size_t LeafSize = 1;
  static const bool PackedTree = true;
  static const bool AllowRemovals = AR;

  bool getUseDivisorCache() const {return true;}
  bool getMinimizeOnInsert() const {return false;}

  bool getDoAutomaticRebuilds() const {return UseDivMask;}
  double getRebuildRatio() const {return 0.5;}
  size_t getRebuildMin() const {return 50;}

private:
  PolyBasis const* mBasis;
  SigPolyBasis const* mSigBasis;
  const Monoid& mMonoid;
  const int mType;
  bool const mPreferSparseReducers;
};

template<class BaseLookup>
class DivLookup;

template<class BL>
class DivLookup : public DivisorLookup {
public:
  typedef BL BaseLookup;
  typedef typename BaseLookup::Monomial Monomial;
  typedef typename BaseLookup::Entry Entry;
  typedef typename BaseLookup::Configuration Configuration;

  static_assert
    (!Configuration::UseTreeDivMask || Configuration::UseDivMask, "");

  DivLookup(const Configuration& c): mLookup(c) {}

  DivLookup(
    const Monoid& monoid,
    int type,
    bool preferSparseReducers
  ):
    mLookup(Configuration(monoid, type, preferSparseReducers))
  {}


  virtual void setBasis(const PolyBasis& basis) {
    mLookup.getConfiguration().setBasis(basis);
  }

  virtual void setSigBasis(const SigPolyBasis& sigBasis) {
    mLookup.getConfiguration().setSigBasis(sigBasis);
  }

  virtual int type() const {return mLookup.getConfiguration().type();}

  template<class Lambda>
  class LambdaWrap {
  public:
    LambdaWrap(Lambda& lambda): mLambda(lambda) {}
    bool proceed(const Entry& entry) const {return mLambda(entry);}
  private:
    Lambda& mLambda;
  };

  template<class Lambda>
  static LambdaWrap<Lambda> lambdaWrap(Lambda& lambda) {
    return LambdaWrap<Lambda>(lambda);
  }

  // *** Signature specific functionality

  virtual size_t regularReducer(ConstMonoRef sig, ConstMonoRef mono) const {
    const auto& conf = mLookup.getConfiguration();
    SigPolyBasis::StoredRatioCmp ratioCmp
      (Monoid::toOld(sig), Monoid::toOld(mono), conf.sigBasis());
    const auto& basis = conf.basis();
    const bool preferSparse = conf.preferSparseReducers();

    auto reducer = size_t(-1);
    auto proceed = [&](const Entry& e) {
      if (ratioCmp.compare(e.index) != GT)
        return true;

      if (reducer != size_t(-1)) {
        if (preferSparse) {
          const auto newTermCount = basis.poly(e.index).nTerms();
          const auto oldTermCount = basis.poly(reducer).nTerms();
          if (newTermCount > oldTermCount)
            return true; // what we already have is sparser
          // resolve ties by picking oldest
          if (newTermCount == oldTermCount && e.index > reducer)
            return true;
        } else { // pick oldest
          if (e.index > reducer)
            return true; // what we already have is older
        }
      }
      reducer = e.index;
      return true;
    };
    auto wrap = lambdaWrap(proceed);
    mLookup.findAllDivisors(mono, wrap);
    return reducer;
  }

  virtual void lowBaseDivisors(
    std::vector<size_t>& divisors,
    size_t maxDivisors,
    size_t newGenerator
  ) const {
    const auto& basis = mLookup.getConfiguration().sigBasis();
    MATHICGB_ASSERT(newGenerator < basis.size());

    auto proceed = [&](const Entry& entry) {
      if (entry.index >= newGenerator)
        return true;
      for (size_t j = 0; j <= divisors.size(); ++j) {
        if (j == divisors.size()) {
          divisors.push_back(entry.index);
          break;
        }
        int cmp = basis.ratioCompare(entry.index, divisors[j]);
        if (cmp == EQ && (entry.index < divisors[j]))
          cmp = GT; // prefer minimum index to ensure deterministic behavior
        if (cmp == GT) {
          divisors.insert(divisors.begin() + j, entry.index);
          break;
        }
      }
      if (divisors.size() > maxDivisors)
        divisors.pop_back();
      MATHICGB_ASSERT(divisors.size() <= maxDivisors);
      return true;
    };
    divisors.clear();
    divisors.reserve(maxDivisors + 1);
    auto wrap = lambdaWrap(proceed);
    mLookup.findAllDivisors(basis.getSignature(newGenerator), wrap);
  }

  virtual size_t highBaseDivisor(size_t newGenerator) const {
    const auto& basis = mLookup.getConfiguration().sigBasis();
    MATHICGB_ASSERT(newGenerator < basis.size());
    auto highDivisor = size_t(-1);

    auto proceed = [&](const Entry& entry) {
      if (entry.index >= newGenerator)
        return true;
      if (highDivisor != size_t(-1)) {
        int cmp = basis.ratioCompare(highDivisor, entry.index);
        if (cmp == LT)
          return true;
        if (cmp == EQ && (entry.index > highDivisor))
          return true; // prefer minimum index to ensure deterministic behavior
      }
      highDivisor = entry.index;
      return true;
    };
    auto wrap = lambdaWrap(proceed);
    mLookup.findAllDivisors(basis.getLeadMonomial(newGenerator), wrap);
    return highDivisor;
  }

  virtual size_t minimalLeadInSig(ConstMonoRef sig) const {
    // Given signature sig, we want to minimize (S/G)g where
    // g and G are the lead term and signature taken over basis elements
    // whose signature G divide S. The code here instead maximizes G/g,
    // which is equivalent and also faster since the basis has a data
    // structure to accelerate comparisons between the ratio of
    // signature to lead term.
    //
    // In case of ties, we select the sparser elements. If there is
    // still a tie, we select the basis element with the largest
    // signature. There can be no further ties since all basis
    // elements have distinct signatures.
    const auto& basis = mLookup.getConfiguration().sigBasis();

    auto minLeadGen = size_t(-1);
    auto proceed = [&](const Entry& entry) {
      if (minLeadGen != size_t(-1)) {
        const int ratioCmp = basis.ratioCompare(entry.index, minLeadGen);
        if (ratioCmp == LT)
          return true;
        if (ratioCmp == EQ) {
          // If same lead monomial in signature, pick the one with fewer terms
          // as that one might be less effort to reduce.
          const size_t minTerms = basis.poly(minLeadGen).nTerms();
          const size_t terms = basis.poly(entry.index).nTerms();
          if (minTerms > terms)
            return true;
          if (minTerms == terms) {
            // If same number of terms, pick the one with larger signature
            // before being multiplied into the same signature. That one
            // might be more reduced as the constraint on regular reduction
            // is less. Also, as no two generators have same signature, this
            // ensures deterministic behavior.
            const auto minSig = basis.getSignature(minLeadGen);
            const auto genSig = basis.getSignature(entry.index);
            const auto sigCmp = basis.monoid().compare(minSig, genSig);
            if (basis.monoid().lessThan(genSig, minSig))
              return true;
          }
        }
      }
      minLeadGen = entry.index;
      return true;
    };
    auto wrap = lambdaWrap(proceed);
    mLookup.findAllDivisors(sig, wrap);
    return minLeadGen;
  }

  // *** Classic GB specific functionality

  virtual size_t classicReducer(ConstMonoRef mono) const {
    const auto& basis = mLookup.getConfiguration().basis();
    const auto preferSparse =
      mLookup.getConfiguration().preferSparseReducers();

    auto reducer = size_t(-1);
    auto proceed = [&](const Entry& entry) {
      if (reducer == size_t(-1)) {
        reducer = entry.index;
        return true;
      }
      if (preferSparse) {
        const auto oldTermCount = basis.poly(reducer).nTerms();
        const auto newTermCount = basis.poly(entry.index).nTerms();
        if (oldTermCount > newTermCount) {
          reducer = entry.index; // prefer sparser
          return true;
        }
        if (oldTermCount < newTermCount)
          return true;
      } // break ties by age

      if (reducer > entry.index)
        reducer = entry.index; // prefer older
      return true;
    };
    auto wrap = lambdaWrap(proceed);
    mLookup.findAllDivisors(mono, wrap);
    return reducer;
  }

  // *** General functionality

  virtual size_t divisor(ConstMonoRef mono) const {
    const Entry* entry = mLookup.findDivisor(mono);
    return entry == 0 ? static_cast<size_t>(-1) : entry->index;
  }

  virtual void divisors(ConstMonoRef mono, EntryOutput& consumer) const {
    auto proceed = [&](const Entry& e) {return consumer.proceed(e.index);};
    auto wrap = lambdaWrap(proceed);
    mLookup.findAllDivisors(mono, wrap);
  }

  virtual void multiples(ConstMonoRef mono, EntryOutput& consumer) const {
    auto proceed = [&](const Entry& e) {return consumer.proceed(e.index);};
    auto wrap = lambdaWrap(proceed);
    mLookup.findAllMultiples(mono, wrap);
  }

  virtual void removeMultiples(ConstMonoRef mono) {
    mLookup.removeMultiples(mono);
  }

  virtual void remove(ConstMonoRef mono) {mLookup.removeElement(mono);}

  virtual size_t size() const {
    return mLookup.size();
  }

  std::string getName() const {return mLookup.getName();}
  const PolyRing& ring() const {return mLookup.configuration().ring();}

  size_t getMemoryUse() const {return mLookup.getMemoryUse();}

  virtual void insert(ConstMonoRef mono, size_t value) {
    mLookup.insert(Entry(mono, value));
  }

private:
  BaseLookup mLookup;
};

MATHICGB_NAMESPACE_END
#endif
