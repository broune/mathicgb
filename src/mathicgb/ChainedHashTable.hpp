// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_CHAINED_HASH_TABLE_GUARD
#define MATHICGB_CHAINED_HASH_TABLE_GUARD

#include <vector>
#include <iostream>
#include <ostream>
#include <memtailor.h>

MATHICGB_NAMESPACE_BEGIN

// One template parameter, with the following:
//   types:
//     ValueType, KeyType
//   functions:
//     size_t HashControl::hash_value(key)
//     bool ::is_equal(key1, key2)
//     void ValueType::combine(value &already_there, value new_one);
//     void ValueType::show(ostream &o, KeyType k, ValueType v)
//
// keys and values: are not copied or stored, except via standard operator=.

class HashControlExample
{
public:
  typedef int * KeyType;
  typedef int ValueType;

  size_t hash_value(KeyType k) const { return static_cast<size_t>(k - static_cast<int *>(0)); }
  bool is_equal(KeyType k1, KeyType k2) const { return k1 == k2; }
  void combine(ValueType &v, ValueType w) const { v += w; }
  void show(std::ostream &o, KeyType k, ValueType v) const { o << "[" << k << " " << v << "]"; }
};

template<typename HashControl>
class ChainedHashTable {
  typedef typename HashControl::KeyType Key;
  typedef typename HashControl::ValueType Value;

  struct node {
    node *next;
    Key key;
    Value value;
    void *unused;
  };

  typedef std::vector<node *> NodeArray;
public:
  struct Stats {
    size_t max_chain_len;
    size_t n_nonempty_bins;
    size_t n_inserts; // # calls to insert a monomial
    size_t n_nodes; // # of unique monomials represented here
    size_t n_eq_true; // total number of true is_equal calls during lookup_and_insert
    size_t n_eq_false;  // same, but number for false.
  };

  ChainedHashTable(HashControl M, int nbits);

  ~ChainedHashTable() {}

  std::string description() const { return "chained hash table"; }

  void reset();  // Clear the table, and memory areas, but don't release space grabbed so far

  void resize(int new_nbits);  // Don't change the nodes, table, but do recreate hashtable_

  bool member(Key k, Value &result_v) const;
  // Return true if k is in the hash table.
  // In this case, set 'result_v' with its value

  void insertUnique(Key m, Value v);
  // The caller must insure that 'm' is not currently in the table.

  bool lookupAndInsert(Key k, Value v);
  // returns true if 'k' is in the table.  In this case, combine values.
  // If false, 'k' is inserted with the given value.

  void getStats(Stats &stats) const; // set the stats table with current values

  void resetStats() const; // reset all stat values to 0

  void dump(int level = 0) const; // For debugging: display the current state of the table

  size_t getMemoryUse() const;
private:
  node * lookup(Key k) const;

  HashControl M_;
  size_t table_size_;
  mutable Stats stats_;

  std::vector<node *> hashtable_;
  memt::Arena mNodes; // where we keep 'node's
};

template <typename HashControl>
ChainedHashTable<HashControl>::ChainedHashTable(HashControl M, int nbits)
  : M_(M),
    table_size_(static_cast<size_t>(1) << nbits)
{
  hashtable_.resize(table_size_);

  // set each entry of hashtable_ to null
  reset();
}

template <typename HashControl>
size_t ChainedHashTable<HashControl>::getMemoryUse() const
{
  size_t total = mNodes.getMemoryUse();
  total += sizeof(node *) * hashtable_.size();
  return total;
}

template <typename HashControl>
void ChainedHashTable<HashControl>::reset()
{
  // Clear the table, and memory areas.
  for (typename NodeArray::iterator i = hashtable_.begin(); i != hashtable_.end(); ++i)
    *i = 0;

  mNodes.freeAllAllocs();
  resetStats();
}

template <typename HashControl>
void ChainedHashTable<HashControl>::resize(int new_nbits)
// Don't change the nodes, table, but do recreate hashtable_
{
  // Make a new vector of node *'s.
  // swap the two.
  // Loop through each one, reinserting the node into the proper bin.

  size_t old_table_size = table_size_;
  table_size_ = static_cast<size_t>(1) << new_nbits;
  NodeArray old_table(table_size_);
  swap(old_table, hashtable_);
  for (size_t i = 0; i<old_table_size; i++)
    {
      node *p = old_table[i];
      while (p != 0)
        {
          node *q = p;
          p = p->next;
          // Reinsert node.  We know that it is unique
          size_t hashval = M_.hash_value(q->key) % table_size_;
          node *r = hashtable_[hashval];
          q->next = r;
          hashtable_[hashval] = q;
        }
    }
}

template <typename HashControl>
bool ChainedHashTable<HashControl>::member(Key m, Value &val) const
{
  node *p = lookup(m);
  if (p == 0) return false;
  val = p->value; // ASSIGN
  return true;
}

template <typename HashControl>
bool ChainedHashTable<HashControl>::lookupAndInsert(Key m, Value val)
// Returns true if m is in the table
{
  node *p = lookup(m);
  if (p == 0)
    {
      insertUnique(m, val);
      return false;
    }
  else
    {
      M_.combine(p->value, val);
      return true;
    }
}

template <typename HashControl>
void ChainedHashTable<HashControl>::insertUnique(Key k, Value v)
// Returns null if not in the table
{
  size_t hashval = M_.hash_value(k) % table_size_; // table_size is a power of 2

  node *q = static_cast<node *>(mNodes.alloc(sizeof(node)));
  q->next = hashtable_[hashval];
  q->key = k; // ASSIGN
  q->value = v; // ASSIGN
  hashtable_[hashval] = q;

  stats_.n_inserts++;
}

template <typename HashControl>
typename ChainedHashTable<HashControl>::node * ChainedHashTable<HashControl>::lookup(Key k) const
// Returns null if not in the table
{
  size_t hashval = M_.hash_value(k) % table_size_; // table_size is a power of 2
  for (node *p=hashtable_[hashval]; p != 0; p = p->next)
    if (M_.is_equal(k, p->key))
      {
        stats_.n_eq_true++;
        return p;
      }
    else
      {
        stats_.n_eq_false++;
      }
  return 0;
}

template <typename HashControl>
void ChainedHashTable<HashControl>::resetStats() const
{
  stats_.max_chain_len = 0;
  stats_.n_nonempty_bins = 0;
  stats_.n_inserts = 0;
  stats_.n_nodes = 0;
  stats_.n_eq_true = 0;
  stats_.n_eq_false = 0;
}

template <typename HashControl>
void ChainedHashTable<HashControl>::getStats(Stats &stats) const
{
  // First we set the values in stats_

  stats_.n_nonempty_bins = 0;
  stats_.max_chain_len = 0;
  stats_.n_nodes = 0;
  for (size_t i = 0; i<table_size_; i++)
    {
      if (hashtable_[i] == 0) continue;
      stats_.n_nonempty_bins++;
      size_t chain_len = 0;
      for (node *p = hashtable_[i]; p != 0; p = p->next)
        chain_len++;
      stats_.n_nodes += chain_len;
      if (chain_len > stats_.max_chain_len)
        stats_.max_chain_len = chain_len;
    }

  if (&stats != &stats_)
    stats = stats_;
}

template <typename HashControl>
void ChainedHashTable<HashControl>::dump(int level) const
{
  // Compute:
  //  # bins in use
  //  max chain length
  //  # keys in use
  //  average number in non-zero bins
  // Report on:
  //  number of is_equal true/false
  //  # of lookup and insert calls
  //
  // For debugging: display the current state of the table
  getStats(stats_);
  std::cout << "ChainedHashTable stats:" << std::endl;
  std::cout << "  # nonempty bins: " << stats_.n_nonempty_bins << std::endl;
  std::cout << "  # nodes: " << stats_.n_nodes << std::endl;
  std::cout << "  max chain length: " << stats_.max_chain_len << std::endl;
  std::cout << "  # insert calls: " << stats_.n_inserts << std::endl;
  std::cout << "  # isequal true calls: " << stats_.n_eq_true << std::endl;
  std::cout << "  # isequal false calls: " << stats_.n_eq_false << std::endl;

  if (level == 0) return;

  for (size_t i = 0; i<table_size_; i++)
    {
      if (hashtable_[i] == 0) continue;
      std::cout << "bin " << i << ": ";
      for (node *p = hashtable_[i]; p != 0; p = p->next)
        M_.show(std::cout, p->key, p->value);
      std::cout << std::endl;
    }
}

MATHICGB_NAMESPACE_END
#endif
