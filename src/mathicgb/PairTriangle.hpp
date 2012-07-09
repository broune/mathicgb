#ifndef _pair_triangle_h
#define _pair_triangle_h

#include "SPairQueue.hpp"
#include "PolyRing.hpp"
class FreeModuleOrder;

typedef unsigned short SmallIndex;
typedef unsigned int BigIndex;

// The following type is only used in the creation of SPairGroups
struct PreSPair {
  BigIndex i;
  monomial signature;
};

class SPairGroup {
public:
  SPairGroup(size_t fixedIndex, monomial signature, std::vector<PreSPair>& prePairs, memt::Arena& arena);

  monomial signature() {return mSignature;}
  void setSignature(monomial signature) {mSignature = signature;}
  const_monomial signature() const {return mSignature;}
  size_t fixedIndex() const {return mFixedIndex;}
  size_t otherIndex() const;

  void increment();
  bool empty() const;
  size_t size() const; // number of S-pairs remaining in this group

  void write(const PolyRing* R, std::ostream& out);

private:
  bool big() const;

  monomial mSignature; // signature of the current pair
  const size_t mFixedIndex; // all pairs here have this index

  union {
    BigIndex* bigBegin; // the current other index is *begin
    SmallIndex* smallBegin;
  };
  union {
    BigIndex* bigEnd; // the other indexes lie in [begin, end)
    SmallIndex* smallEnd;
  };
};


// Object that stores integer pairs.
// The pairs are stored in the shape of a triangle:
//
//          3
//        2 2
//      1 1 1
//    0 0 0 0
//  ---------
//  0 1 2 3 4
//
// So column a stores the numbers b such that (a,b) is a pair
// with a > b. Not all pairs need be present, and there can be
// any ordering of the entries in each column.
//
// PairTriangles are useful for storing S-pair indices, and the
// class has code that helps with that.
class PairTriangle {
public:
  PairTriangle(
    const FreeModuleOrder& order,
    const PolyRing& ring,
    size_t queueType);
  ~PairTriangle();

  // Returns how many columns the triangle has
  size_t columnCount() const {return mColumnCount;}

  // Returns how many pairs are in the triangle
  size_t size() const;

  // Returns true if there are no pairs in the triangle
  bool empty() const {return mQueue->empty();}

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

protected:
  // Sub classes implement this to say what monomial each pair is ordered
  // according to. That monomial should be placed into orderBy.
  //
  // If false is returned, the requested S-pair is not valid and should be
  // skipped.
  virtual bool calculateOrderBy(size_t a, size_t b, monomial orderBy) const = 0;

private:
  bool const mUseSingletonGroups;
  size_t mColumnCount;
  std::vector<SPairGroup*> mGroups;
  std::auto_ptr<SPairQueue> mQueue;
  memt::Arena mArena;
  FreeModuleOrder const& mOrder;
  std::vector<PreSPair> mPrePairs;
  std::vector<PreSPair> mPrePairTmp;
  PolyRing const& mRing;
};

#endif
