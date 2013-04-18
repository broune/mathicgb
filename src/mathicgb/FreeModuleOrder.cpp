// Copyright 2011 Michael E. Stillman

#include "stdinc.h"
#include "FreeModuleOrder.hpp"

#include "Poly.hpp"
#include "Ideal.hpp"
#include "SigSPairQueue.hpp"
#include "GroebnerBasis.hpp"
#include "PolyRing.hpp"
#include <mathic.h>
#include <iostream>
#include <algorithm>
#include <limits>
#include <stdexcept>

// In this file we define various term orders. Some classes need to be
// specialized to each term order which we do using templates. Those
// classes are defined here since they should not be used from
// anywhere else - instead they should be accessed through their
// virtual interface.

// *******************************************************
// ** Utility objects

namespace {
  // Comparer for scrambled signatures
  template<class Cmp>
  class ScrambledComparer {
  public:
    ScrambledComparer(const Cmp& cmp): mCmp(cmp), mComparisons(0) {}
    bool operator()(const PreSPair& a, const PreSPair& b) const {
      ++mComparisons;
      return mCmp.scrambledLessThan(a.signature, b.signature);
    }
    size_t comparisons() const {return mComparisons;}
  private:
    const Cmp& mCmp;
    mutable size_t mComparisons;
  };

  // Iterator that accesses the field i based on a passed-in iterator.
  template<class PairIterator>
  class IndexIterator {
  public:

	typedef typename PairIterator::iterator_category iterator_category;
    typedef decltype(reinterpret_cast<typename PairIterator::value_type*>(0)->i) value_type;
	typedef typename PairIterator::difference_type difference_type;
	typedef value_type* pointer;
	typedef value_type& reference;

	IndexIterator(PairIterator pairIterator): mIterator(pairIterator) {}
	IndexIterator& operator++() {++mIterator; return *this;}
    const value_type operator*() const {return mIterator->i;}
	difference_type operator-(const IndexIterator<PairIterator>& it) const {
	  return mIterator - it.mIterator;
	}
	bool operator==(const IndexIterator<PairIterator>& it) const {
	  return mIterator == it.mIterator;
	}
	bool operator!=(const IndexIterator<PairIterator>& it) const {
	  return mIterator != it.mIterator;
	}

  private:
	PairIterator mIterator;
  };
}

// *******************************************************
// ** SigSPairQueue
namespace {
  // Configuration of mathic::PairTriangle for use with signature queues.
  template<class Cmp>
  class SigPTConfiguration {
  public:
	SigPTConfiguration
    (GroebnerBasis const& basis, Cmp const& cmp): mBasis(basis), mCmp(cmp) {}

	typedef monomial PairData;
	void computePairData(size_t col, size_t row, monomial sig) {
      MATHICGB_ASSERT(mBasis.ratioCompare(col, row) != EQ);
      // ensure that ratio(col) > ratio(row)
      if (mBasis.ratioCompare(col, row) == LT)
        std::swap(col, row);
      mBasis.ring().monomialFindSignature
        (mBasis.getLeadMonomial(col),
         mBasis.getLeadMonomial(row),
         mBasis.getSignature(col), sig);
      mCmp.scrambleSignatureForFastComparison(sig);
    }

	typedef bool CompareResult;
	bool compare(int colA, int rowA, const_monomial a,
				 int colB, int rowB, const_monomial b) const {
      return mCmp.scrambledLessThan(b, a);
	}
	bool cmpLessThan(bool v) const {return v;}

    GroebnerBasis const& basis() const {return mBasis;}
    Cmp const& comparer() const {return mCmp;}

  private:
    GroebnerBasis const& mBasis;
    Cmp const& mCmp;
  };

  // Object that stores S-pairs and orders them according to a monomial
  // or signature.
  template<class Cmp>
  class ConcreteSigSPairQueue : public SigSPairQueue {
  public:
    ConcreteSigSPairQueue(GroebnerBasis const& basis, Cmp const& cmp):
      mPairQueue(PC(basis, cmp)) {}

    virtual monomial popSignature(Pairs& pairs) {
      pairs.clear();
      if (mPairQueue.empty())
        return 0;
      monomial sig = ring().allocMonomial();
      ring().monomialCopy(mPairQueue.topPairData(), sig);
      do {
        pairs.push_back(mPairQueue.topPair());
        mPairQueue.pop();
      } while
          (!mPairQueue.empty() && ring().monomialEQ(mPairQueue.topPairData(), sig));
      comparer().unscrambleSignature(sig);
      return sig;
    }

    virtual void pushPairs(size_t pairWith, IndexSigs& pairs) {
#ifdef DEBUG
      monomial tmp = ring().allocMonomial();
      for (size_t i = 0; i < pairs.size(); ++i) {
        MATHICGB_ASSERT(pairs[i].i < columnCount());
        mPairQueue.configuration().computePairData
          (columnCount(), pairs[i].i, tmp);
        comparer().unscrambleSignature(tmp);
        MATHICGB_ASSERT(ring().monomialEQ(tmp, pairs[i].signature));
      }
      ring().freeMonomial(tmp);
#endif

      if (columnCount() >= std::numeric_limits<BigIndex>::max())
        throw std::overflow_error
          ("Too large basis element index in constructing S-pairs.");
    
      // sort and insert new column
      ScrambledComparer<Cmp> cmp(comparer());
      comparer().scrambleSignaturesForFastComparison(pairs);
      std::sort(pairs.begin(), pairs.end(), cmp);
      //mPreComparisons += cmp.comparisons();
      //order().destructiveSort(pairs);
      typedef IndexIterator<std::vector<PreSPair>::const_iterator> Iter;
      mPairQueue.addColumnDescending(Iter(pairs.begin()), Iter(pairs.end()));
    
      // free signatures
      std::vector<PreSPair>::iterator end = pairs.end();
      for (std::vector<PreSPair>::iterator it = pairs.begin(); it != end; ++it)
        ring().freeMonomial(it->signature);
      pairs.clear();
    }
  
    virtual std::string name() const {return "todo";}

    virtual size_t memoryUse() const {return mPairQueue.getMemoryUse();}

    virtual size_t pairCount() const {return mPairQueue.pairCount();}

    virtual size_t columnCount() const {return mPairQueue.columnCount();}

  private:
    ConcreteSigSPairQueue(const ConcreteSigSPairQueue<Cmp>&); // not available
    void operator=(const ConcreteSigSPairQueue<Cmp>&); // not available

    // the compiler should be able to resolve these accessors into a direct
    // offset as though these were member variables.
    const GroebnerBasis& basis() const {
      return mPairQueue.configuration().basis();
    }
    const Cmp& comparer() const {return mPairQueue.configuration().comparer();}
    PolyRing const& ring() const {return basis().ring();}
    FreeModuleOrder const& order() const {return basis().order();}

    typedef SigPTConfiguration<Cmp> PC;
    mathic::PairQueue<PC> mPairQueue;
    friend struct mathic::PairQueueNamespace::ConstructPairDataFunction<PC>;
    friend struct mathic::PairQueueNamespace::DestructPairDataFunction<PC>;
  };
}

namespace mathic {
  namespace PairQueueNamespace {
	template<class Cmp>
    struct ConstructPairDataFunction<SigPTConfiguration<Cmp> > {
      inline static void function
      (void* memory, Index col, Index row, SigPTConfiguration<Cmp>& conf) {
        MATHICGB_ASSERT(memory != 0);
        MATHICGB_ASSERT(col > row);
        monomial* pd = new (memory) monomial
          (conf.basis().ring().allocMonomial());
        conf.computePairData(col, row, *pd);
      }
	};

	template<class Cmp>
    struct DestructPairDataFunction<SigPTConfiguration<Cmp> > {
      inline static void function
      (monomial* pd, Index col, Index row, SigPTConfiguration<Cmp>& conf) {
        MATHICGB_ASSERT(pd != 0);
        MATHICGB_ASSERT(col > row);
        conf.basis().ring().freeMonomial(*pd);
      }
    };
  }
}

// *******************************************************
// ** Term orders

template<class Cmp>
class ConcreteOrder : public FreeModuleOrder {
public:
  ConcreteOrder(Cmp cmp): mCmp(cmp), mComparisons(0), mPreComparisons(0) {}
  virtual ~ConcreteOrder() {}

  virtual int signatureCompare(const_monomial sigA, const_monomial sigB) const {
    ++mComparisons;
    return mCmp.signatureCompare(sigA, sigB);
  }

  virtual int signatureCompare(
    const_monomial sigA,
    const_monomial monoB,
    const_monomial sigB
  ) const {
    ++mComparisons;
    return mCmp.signatureCompare(sigA, monoB, sigB);
  }

  virtual void sortAndScrambleSignatures(std::vector<PreSPair>& pairs) const {
    ScrambledComparer<Cmp> cmp(mCmp);
    mCmp.scrambleSignaturesForFastComparison(pairs);
    std::sort(pairs.begin(), pairs.end(), cmp);
    mPreComparisons += cmp.comparisons();
  }

  virtual void appendBasisElement(const_monomial m) {
    mCmp.appendBasisElement(m);
  }

  virtual std::string description() const {
    return mCmp.description();
  }

  virtual void getStats(size_t& comparisons, size_t& preComparisons) const {
    comparisons = mComparisons;
    preComparisons = mPreComparisons;
  }

  virtual std::unique_ptr<SigSPairQueue>
  createSigSPairQueue(GroebnerBasis const& basis) const {
    return std::unique_ptr<SigSPairQueue>
      (new ConcreteSigSPairQueue<Cmp>(basis, mCmp));
  }

private:
  Cmp mCmp;
  mutable size_t mComparisons;
  mutable size_t mPreComparisons;
};



// ** Graded reverse lex.
// Degrees and exponents considered from high index to low index.
//
//   rule: a[i] < b[i] then a > b
//
// also applies to component, which is considered last.
class OrderA {
public:
  OrderA(const PolyRing* ring): mRing(ring) {}

  int signatureCompare(const_monomial sigA, const_monomial sigB) const {
    return mRing->monomialCompare(sigA, sigB);
  }

  void scrambleSignatureForFastComparison(monomial sig) const {}
  void scrambleSignaturesForFastComparison(std::vector<PreSPair>& pairs) const {}
  bool scrambledLessThan(const_monomial sigA, const_monomial sigB) const {
    return mRing->monomialLT(sigA, sigB);
  }
  void unscrambleSignature(monomial sig) const {}

  int signatureCompare(const_monomial sigA,
    const_monomial monoB, const_monomial sigB) const {
    return mRing->monomialCompare(sigA, monoB, sigB);
  }

  void appendBasisElement(const_monomial) {}

  char const* description() const {return "GrevLex IndexDown";}

private:
  PolyRing const* const mRing;
};

class OrderC
{
public:
  OrderC(const Ideal* I):
    mRing(I->getPolyRing()),
    topindex(I->getPolyRing()->maxMonomialSize() - 2)
  {}

  void appendBasisElement(const_monomial m) {}

  void scrambleSignatureForFastComparison(monomial sig) const {}
  void scrambleSignaturesForFastComparison(std::vector<PreSPair>& pairs) const {}
  bool scrambledLessThan(const_monomial a, const_monomial b) const {
    return signatureCompare(a, b) == LT;
  }
  void unscrambleSignature(monomial sig) const {}

  int signatureCompare(const_monomial sig, const_monomial sig2) const {
    int da = - sig[topindex];
    int db = -sig2[topindex];
    if (da > db) return GT;
    if (db > da) return LT;
    int cmp = *sig  - *sig2;
    if (cmp < 0) return GT;
    if (cmp > 0) return LT;

    auto a = sig;
    auto b = sig2;
    for (size_t i = topindex; i >= 1; i--) {
      int cmp = a[i] - b[i];
      if (cmp != 0)
        return cmp < 0 ? GT : LT;
    }
    return EQ;
  }

  int signatureCompare(
    const_monomial sig,
    const_monomial m2,
    const_monomial sig2
  ) const {
    int da = - sig[topindex];
    int db = - m2[topindex];
    db += -sig2[topindex];
    if (da > db) return GT;
    if (db > da) return LT;
    int cmp = *sig  - *sig2;
    if (cmp < 0) return GT;
    if (cmp > 0) return LT;

    auto a = sig;
    auto b = sig2;
    for (size_t i = topindex; i >= 1; i--)
      {
        int cmp = a[i] - b[i] - m2[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
    return EQ;
  }

  char const* description() const {return "DegreeUp IndexDown GrevLex";}

private:
  PolyRing const* const mRing;
  const size_t topindex;
};

// Let l(ae_i) be the leading monomial of ag_i.
// Indexes are considered from high to low
// rule 1: higher component wins (if mUp is true)
// rule 1': lower components wins (if mUp is false)
// rule 2: higher degree of l(ae_I) wins
// rule 3: reverse lex on l(ae_i)
class OrderE
{
public:
  OrderE(Ideal const* I):
    mRing(I->getPolyRing()),
    topindex(mRing->maxMonomialSize() - 2)
  {}

  void appendBasisElement(const_monomial m) {}

  int signatureCompare(const_monomial a, const_monomial b) const {
    if (*a != *b)
      return *a < *b ? GT : LT;
    for (size_t i = topindex; i >= 1; i--) {
      int cmp = a[i] - b[i];
      if (cmp != 0)
        return cmp < 0 ? GT : LT;
    }
    return EQ;
  }

  void scrambleSignatureForFastComparison(monomial sig) const {}
  void scrambleSignaturesForFastComparison(std::vector<PreSPair>& pairs) const {}
  inline bool scrambledLessThan(const_monomial a, const_monomial b) const {
    return signatureCompare(a,b) == LT;
  }
  void unscrambleSignature(monomial sig) const {}

  int signatureCompare(const_monomial a, const_monomial m2, const_monomial b) const {
    int cmp = *a - *b;
    if (cmp != 0) {
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }

    for (size_t i = topindex; i >= 1; i--) {
      int cmp = a[i] - b[i] - m2[i];
      if (cmp < 0) return GT;
      if (cmp > 0) return LT;
    }
    return EQ;
  }

  virtual char const* description() const {
    return "IndexDown SchreyerGrevLex";
  }

private:
  PolyRing const* const mRing;
  size_t const topindex; // taken from R
};

void FreeModuleOrder::displayOrderTypes(std::ostream &o)
{
  o << "FreeModule orders:" << std::endl;
  o << "  1   GrevLex IndexDown" << std::endl;
  o << "  2   DegreeUp IndexUp GrevLex" << std::endl; // done
  o << "  3   DegreeUp IndexDown GrevLex" << std::endl;
  o << "  4   SchreyerGrevLexUp" << std::endl; // done
  o << "  5   SchreyerGrevLexDown" << std::endl; // done
  o << "  6   IndexUp SchreyerGrevLex" << std::endl; // done
  o << "  7   IndexDown SchreyerGrevLex" << std::endl;
}

std::unique_ptr<FreeModuleOrder> FreeModuleOrder::makeOrder(FreeModuleOrderType type, const Ideal* I)
{
  if (type == 0)
    type = 1;  // Set the default

  switch (type) {
  case 4:
  case 5:
  case 1:
    return make_unique<ConcreteOrder<OrderA>>(OrderA(I->getPolyRing()));

  case 2:
  case 3:
   return make_unique<ConcreteOrder<OrderC>>(OrderC(I));

  case 6:
  case 7:
    return make_unique<ConcreteOrder<OrderE>>(OrderE(I));

  default: break;
  }

  std::cerr << "unknown free module order type" << std::endl;
  std::cerr << "possible orders are: " << std::endl;
  for (size_t i = 1; i <= 7; ++i) {
    auto order = makeOrder(static_cast<FreeModuleOrderType>(i), I);
    std::cerr << "  " << i << ": " << order->description() << std::endl;
  }
  exit(1);
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
