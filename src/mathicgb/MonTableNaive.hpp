// Copyright 2011 Michael E. Stillman

#ifndef _monomial_table_h_
#define _monomial_table_h_

#include <memtailor.h>
#include "PolyRing.hpp"

struct mon_node { // each node is in 'nodes' arena
  mon_node *next;
  const_monomial monom; // points into 'monoms'
};

class MonTableNaiveConfiguration
{
public:
  MonTableNaiveConfiguration(const PolyRing *R) : R_(R) {}

  const PolyRing * getPolyRing() const { return R_; }

private:
  const PolyRing *R_;
};

class MonTableNaive {
public:
  typedef size_t ValueType;
  typedef MonTableNaiveConfiguration Configuration;

  MonTableNaive(const PolyRing *R);
  MonTableNaive(const Configuration& C);
  ~MonTableNaive();

  Configuration& getConfiguration() { return conf_; }
  const Configuration& getConfiguration() const { return conf_; }

  const PolyRing * getPolyRing() const { return conf_.getPolyRing(); }

  bool member(const_monomial t, ValueType & val) const;
  bool insert(const_monomial t, ValueType val = 0); // Only insert if not there

  std::string getName() const;
  void display(std::ostream &o, int level) const;
  void dump(int level) const;

  size_t n_elems() const;

  void displayStats(std::ostream &o) const;

  struct Stats {
    size_t n_member;
    size_t n_inserts;  // includes koszuls
    size_t n_insert_already_there;
    size_t n_compares;
  };

  void getStats(Stats &stats) const { stats = stats_; }

  size_t getMemoryUse() const;

  void getMonomials(std::vector<const_monomial>& monomials);

private:
  void insert_node(mon_node *p, const_monomial t);

  // For now, keep a sorted linked list of monomials (in increasing order), one such list for each component found
  // member: just goes through each one, checking for divisibility
  // insert: find insertion spot, insert it, then try to remove ones after that that are not needed
  Configuration conf_;

  memt::BufferPool mPool; // used for mon_node's
  mon_node *table;

  // stats:
  mutable Stats stats_;
};


#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
