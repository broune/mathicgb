#ifndef _pair_triangle_h
#define _pair_triangle_h

#include <memtailor.h>
#include <mathic.h>
#include "PolyRing.hpp"
#include "FreeModuleOrder.hpp"
class FreeModuleOrder;

/*typedef unsigned short SmallIndex;
typedef unsigned int BigIndex;

// The following type is only used in the creation of SPairGroups
struct PreSPair {
  BigIndex i;
  monomial signature;
  };*/

// Object that stores S-pairs and orders them according to a monomial
// or signature.
class PairTriangle {
public:
  PairTriangle(
    const FreeModuleOrder& order,
    const PolyRing& ring,
    size_t queueType);

  // Returns how many columns the triangle has
  size_t columnCount() const {return mPairQueue.columnCount();}

  // Returns how many pairs are in the triangle
  size_t pairCount() const {return mPairQueue.pairCount();}

  // Returns true if there are no pairs in the triangle
  bool empty() const {return mPairQueue.empty();}

  // Adds a new column of the triangle and opens it for addition of pairs.
  // This increases columnCount() by one, and the index of the new column
  // is the previous value of columnCount(). Must call endColumn
  // before calling beginColumn again or using the new column.
  void beginColumn();

  // Adds a pair to the most recent column that must still be open for
  // addition of pairs. If a is the index of the new column, then
  // the added pair is (a,index). index must be less than a.
  // orderBy must have been allocated on the ring's pool of monomials,
  // and ownership of the memory is passed to the this triangle object.
  void addPair(size_t index, monomial orderBy);

  // Closes the new column for addition of pairs. Must be preceded by a call
  // to beginColumn(). This sorts the added pairs according to their orderBy
  // monomials.
  void endColumn();

  // Returns a pair with minimal orderBy monomial.
  std::pair<size_t, size_t> topPair() const;

  // Returns the minimal orderBy monomial along all pairs. This is the orderBy
  // monomial of topPair().
  const_monomial topOrderBy() const;

  // Removes topPair() from the triangle.
  void pop();

  size_t getMemoryUse() const;

  std::string name() const;
  size_t mColumnCount;

protected:
  // Sub classes implement this to say what monomial each pair is ordered
  // according to. That monomial should be placed into orderBy.
  //
  // If false is returned, the requested S-pair is not valid and should be
  // skipped.
  virtual bool calculateOrderBy(size_t a, size_t b, monomial orderBy) const = 0;

private:
  FreeModuleOrder const& mOrder;
  std::vector<PreSPair> mPrePairs;
  PolyRing const& mRing;

  class PC {
  public:
	PC(PairTriangle const& tri): mTri(tri) {}
	
	typedef monomial PairData;
	void computePairData(size_t col, size_t row, monomial m) {
	  mTri.calculateOrderBy(col, row, m);
	}

	typedef bool CompareResult;
	bool compare(int colA, int rowA, const_monomial a,
				 int colB, int rowB, const_monomial b) const {
	  return mTri.mOrder.signatureCompare(a, b) == GT;
	}
	bool cmpLessThan(bool v) const {return v;}

	// these are not required for a configuration but we will use
	// them from this code.
	monomial allocPairData() {return mTri.mRing.allocMonomial();}
	void freePairData(monomial m) {mTri.mRing.freeMonomial(m);}

  private:
	PairTriangle const& mTri;
  };
  mathic::PairQueue<PC> mPairQueue;
  friend void mathic::PairQueueNamespace::constructPairData<PC>(void*,Index,Index,PC&);
  friend void mathic::PairQueueNamespace::destructPairData<PC>(monomial*,Index,Index, PC&);
};

namespace mathic {
  namespace PairQueueNamespace {
	template<>
	inline void constructPairData<PairTriangle::PC>
	(void* memory, Index col, Index row, PairTriangle::PC& conf) {
	  MATHICGB_ASSERT(memory != 0);
	  MATHICGB_ASSERT(col > row);
	  monomial* pd = new (memory) monomial(conf.allocPairData());
	  conf.computePairData(col, row, *pd);
	}

	template<>
	inline void destructPairData
	(monomial* pd, Index col, Index row, PairTriangle::PC& conf) {
	  MATHICGB_ASSERT(pd != 0);
	  MATHICGB_ASSERT(col > row);
	  conf.freePairData(*pd);
	}	
  }
}

#endif
