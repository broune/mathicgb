// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_POLY_HASH_TABLE_GUARD
#define MATHICGB_POLY_HASH_TABLE_GUARD

#include "PolyRing.hpp"
#include "Poly.hpp"
#include <utility>
#include <memtailor.h>
#include <vector>

MATHICGB_NAMESPACE_BEGIN

// The hash table is a map:  monomial => coeff
// Operations required on monomials:
//  hash (this will currently pick out one entry of a monomial)
//
// Operations required on coefficients:
//  add, maybe multiply too?
//  isZero?
//
// Does not take ownership of any of the monomials.
class PolyHashTable {
public:
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::MonoRef MonoRef;
  typedef Monoid::ConstMonoRef ConstMonoRef;
  typedef Monoid::MonoPtr MonoPtr;
  typedef Monoid::ConstMonoPtr ConstMonoPtr;
  typedef coefficient Value;

  class Node {
  public:

    const_monomial& mono() {return monom;}
    const const_monomial& mono() const {return monom;}

    Value& value() {return coeff;}
    const Value& value() const {return coeff;}

  private:
    friend class PolyHashTable;

    Node*& next() {return mNext;}
    Node* next() const {return mNext;}

    Node* mNext;
    coefficient coeff;
    const_monomial monom;
  };
  typedef Node node; // todo: remove

  typedef std::vector<node*> MonomialArray;

  PolyHashTable(const PolyRing *R, int nbits);

  std::string description() const {return "polynomial hash table";}

  void reset();  // Clear the table, and memory areas.

  void resize(size_t new_nbits);  // Don't change the nodes, table, but do recreate hashtable_

  size_t getMemoryUse() const;

  //@ insert multiplier * g: any monomials already in the hash table are removed,
  // but their field coefficients in the table are modified accordingly.
  // Resulting pointers to 'node's are placed, in order, into result.
  void insert(Poly::const_iterator first, 
              Poly::const_iterator last,
              MonomialArray &result);
  void insert(const_term multiplier, 
              Poly::const_iterator first, 
              Poly::const_iterator last,
              MonomialArray &result);
  void insert(const_monomial multiplier, 
              Poly::const_iterator first, 
              Poly::const_iterator last,
              MonomialArray &result);

  // Inserts t into the hashtable. Returns true if there is not already
  // a node with the monomial t.monom. If there is already such a node,
  // then t.coeff is added to the coefficient of that node. In either case,
  // node will point to the node for t.monom. The original value of
  // nodeOut is not used.
  std::pair<bool, node*> insert(const_term termToInsert);

  // Removes the node from the hash table.
  void remove(node* n);

  // deprecated: use remove instead
  // popTerm removes 'p' from the hash table, setting result_coeff and result_monom if the coefficient is not zero.
  // result_monom is set with a pointer into monomial space in this class, so will only
  // be valid until a 'reset' is called.
  bool popTerm(node *p, coefficient &result_coeff, const_monomial &result_monom);

protected:
  node * makeNode(coefficient coeff, const_monomial monom);
  void unlink(node *p);
  bool lookup_and_insert(const_monomial m, coefficient val, node *&result);

  const PolyRing& mRing;
  std::vector<node*> mHashTable;
  size_t mHashMask; // this is the number, in binary:  00001111...1, where
                    // the number of 1's is mLogTableSize

  memt::Arena mArena; // space for monomials represented in this class.  Also nodes??

  size_t mTableSize;
  size_t mLogTableSize; // number of bits in the table: mTableSize should be 2^mLogTableSize

  size_t mNodeCount;  // current number of nodes in the hash table
  size_t mBinCount;
  size_t mMaxCountBeforeRebuild;

  size_t mMonomialSize;
};

MATHICGB_NAMESPACE_END
#endif
