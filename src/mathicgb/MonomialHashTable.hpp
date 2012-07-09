// Copyright 2011 Michael E. Stillman

#ifndef _MonomialHashTable_h_
#define _MonomialHashTable_h_

#include "ChainedHashTable.hpp"
#include "PolyRing.hpp"

class MonomialHashControl
{
public:
  typedef const_monomial KeyType;
  typedef int ValueType;

  MonomialHashControl(const PolyRing *R) : R_(R), hash_index_(R->monomialHashIndex()) {}
  size_t hash_value(KeyType k) const { return R_->monomialHashValue(k); }
  bool is_equal(KeyType k1, KeyType k2) const { return R_->monomialEQ(k1, k2); }
  void combine(ValueType &, ValueType) const { }
  void show(std::ostream &o, KeyType k, ValueType v) const { 
    o << "["; 
    R_->monomialDisplay(o, k);
    o << " " << v << "]"; 
  }
private:
  const PolyRing *R_;
  size_t hash_index_;
};

typedef ChainedHashTable<MonomialHashControl> MonomialHashTableBasic;

class MonomialHashTable
{
public:
  typedef MonomialHashTableBasic::Stats Stats;

  MonomialHashTable(const PolyRing *R, int nbits)
  : R_(R), H_(R, nbits) {}

  ~MonomialHashTable() {}

  std::string description() const { return "monomial hash table"; }

  void reset() {
    // Clear the table, and memory areas, but don't release space grabbed so far
    mMonomialPool.freeAllAllocs();
    H_.reset();
  }

  void resize(int new_nbits) { H_.resize(new_nbits); }
  // Don't change the nodes, table, but do recreate hashtable_

  bool member(const_monomial m, int &result_val) const { return H_.member(m, result_val); }
    // Return true if m is in the hash table.
    // In this case, set 'result_val' with its value

  void insertUnique(const_monomial m, int v) {
    // The caller must insure that 'm' is not currently in the table.
    // First copy m.
    monomial m1 = R_->allocMonomial(mMonomialPool);
    R_->monomialCopy(m, m1);
    H_.insertUnique(m1,v);
  }

  bool lookupAndInsert(const_monomial m, int  &v) {
    // returns true if 'm' is in the table.  In this case, v is set to value that is in the table,
    // If false, 'm' is inserted, and its value is set to v.
    if (member(m, v)) return true;
    insertUnique(m, v);
    return false;
  }

  void getStats(Stats &stats) const { H_.getStats(stats) ; }

  void resetStats() const { H_.resetStats(); }
  // reset all values to 0

  void dump(int level = 0) const { H_.dump(level); }  // For debugging: display the current state of the table

  size_t getMemoryUse() const { return mMonomialPool.getMemoryUse() + H_.getMemoryUse(); }
private:
  const PolyRing *R_;
  memt::Arena mMonomialPool;
  MonomialHashTableBasic H_;
};
#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
