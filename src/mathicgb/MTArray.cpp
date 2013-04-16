#include "stdinc.h"
#include "MTArray.hpp"

#include "MonTableNaive.hpp"
#include "MonTableKDTree.hpp"
#include "MonTableDivList.hpp"

template <typename MT>
class MTArrayT : public MonomialTableArray
{
  typedef MT T;
  typedef typename MT::Configuration Conf;
  typedef MonomialTableArray::V V; // value type
public:
  MTArrayT(size_t components, const Conf &conf):
    MonomialTableArray(conf.getPolyRing()), conf_(conf) {
    for (size_t i = 0; i < components; ++i)
      addComponent();
  }

  ~MTArrayT()
  {
    for (size_t i = 0; i<tables.size(); i++)
      delete tables[i];
  }

  Conf& getConfiguration() { return conf_; }
  const Conf &getConfiguration() const { return conf_; }

  bool insert(const_monomial m, V val);
  // Only inserts if minimal, in this case "true" is returned.
  // and the object takes ownership of 'm'.

  void addComponent() {tables.push_back(new T(conf_));}

  bool member(const_monomial m, V &result);

  void getStats(MonomialTableArray::Stats &stats) const;

  std::string description() const;

  void displayStats(std::ostream &o) const;

  void display(std::ostream &o, int level) const;

  virtual void printFrobbyM2Format
    (std::ostream& out, size_t component) const;

  void dump(int level=0) const;

  size_t n_elems() const;

  size_t getMemoryUse() const;
private:
  Conf conf_; // Used to create new instances of T.
  std::vector<T *> tables;

  mutable MonomialTableArray::Stats stats_;
};

template <typename MT>
bool MTArrayT<MT>::insert(const_monomial m, V val)
{
  stats_.n_calls_insert++;
  size_t x = R->monomialGetComponent(m);
  MATHICGB_ASSERT(x < tables.size());
  MATHICGB_ASSERT(tables[x] != 0);
  bool result = tables[x]->insert(m, val);
  if (result) stats_.n_actual_inserts++;
  return result;
}

template <typename MT>
bool MTArrayT<MT>::member(const_monomial m, V &result)
{
  stats_.n_calls_member++;
  size_t x = R->monomialGetComponent(m);
  MATHICGB_ASSERT(x < tables.size());
  MATHICGB_ASSERT(tables[x] != 0);
  return tables[x]->member(m, result);
}

template <typename MT>
void MTArrayT<MT>::getStats(Stats &stats) const
{
  stats = stats_;

  // compute denseness
  std::vector<const_monomial> monomials;
  unsigned long long support = 0;
  unsigned long long exponentCount = 0;
  for (size_t i=0; i<tables.size(); i++) {
    T *p = tables[i];
    MATHICGB_ASSERT(p != 0);
    monomials.clear();
    p->getMonomials(monomials);
    typedef std::vector<const_monomial>::const_iterator iter;
    for (iter it = monomials.begin(); it != monomials.end(); ++it) {
      support += R->monomialSizeOfSupport(*it);
      exponentCount += R->getNumVars();
    }
  }
  stats.denseness = static_cast<double>(support) / exponentCount;
}

template <typename MT>
void MTArrayT<MT>::displayStats(std::ostream & /* o */) const
{
}

template <typename MT>
void MTArrayT<MT>::display(std::ostream &o, int level) const
{
  std::vector<const_monomial> monomials;
  for (size_t i=0; i<tables.size(); i++)
    {
      T *p = tables[i];
      MATHICGB_ASSERT(p != 0);
      if (p->n_elems() == 0)
        continue;

      o << "  " << i << ": ";
      monomials.clear();
      p->getMonomials(monomials);
      std::sort(monomials.begin(), monomials.end(),
        MonomialTableArray::MonomialCompare(*R));
      typedef std::vector<const_monomial>::const_iterator iter;
      for (iter it = monomials.begin(); it != monomials.end(); ++it)
        {
          R->monomialDisplay(o, *it, false);
          o << "  ";
        }
      o << '\n';
    }
}

template <typename MT>
void MTArrayT<MT>::printFrobbyM2Format
  (std::ostream& out, size_t component) const
{
  MATHICGB_ASSERT(component < tables.size());
  T* p = tables[component];
  MATHICGB_ASSERT(p != 0);
  std::vector<const_monomial> monomials;
  p->getMonomials(monomials);
  std::sort(monomials.begin(), monomials.end(),
    MonomialTableArray::MonomialCompare(*R));

  out << "I = monomialIdeal(\n";
  typedef std::vector<const_monomial>::const_iterator iter;
  for (iter it = monomials.begin(); it != monomials.end(); ++it) {
    if (it != monomials.begin())
      out << ",\n";
    R->printMonomialFrobbyM2Format(out, *it);
  }
  out << "\n);\n";
}

template <typename MT>
void MTArrayT<MT>::dump(int level) const
{
  // display on stderr the table.
  displayStats(std::cerr);
  if (level > 0) display(std::cerr, level-1);
}

template <typename MT>
size_t MTArrayT<MT>::n_elems() const
{
  size_t result = 0;
  for (size_t i=0; i<tables.size(); i++)
    if (tables[i] != 0)
      result += tables[i]->n_elems();
  return result;;
}

template <typename MT>
size_t MTArrayT<MT>::getMemoryUse() const
{
  size_t sum = tables.capacity() * sizeof(tables.front());
  for (size_t i = 0; i < tables.size(); ++i) {
    MATHICGB_ASSERT(tables[i] != 0);
    sum += sizeof(tables[i]) + tables[i]->getMemoryUse();
  }
  return sum;
}

template <typename MT>
std::string MTArrayT<MT>::description() const
{
  T obj(getConfiguration());
  return obj.getName();
}

int MonomialTableArray::displayMTTypes(std::ostream &o)
 // returns n s.t. 0..n-1 are valid types
{
  o << "Monomial table types:" << std::endl;
  o << "  0   naive" << std::endl;
  o << "  1   divlist+divmask" << std::endl;
  o << "  2   kdtree+divmask" << std::endl;
  o << "  3   divlist" << std::endl;
  o << "  4   kdtree" << std::endl;
  return 5;
}

std::unique_ptr<MonomialTableArray> MonomialTableArray::make(const PolyRing *R, int typ, size_t components, bool allowRemovals)
{
  switch (typ) {
    case 0:
      return make_unique<MTArrayT<MonTableNaive>>(components, R);

    case 1: {
      typedef MonTableDivList<false,true> MT;
      return make_unique<MTArrayT<MT>>(components, MT::Configuration
        (R, true, false, true, .5, 500));
    }
    case 2: {
      if (allowRemovals) {
        typedef MonTableKDTree<true,true,1,true> MT;
        return make_unique<MTArrayT<MT>>(components, MT::Configuration
          (R, true, false, true, .5, 50));
      } else {
        typedef MonTableKDTree<true,true,1,false> MT;
        return make_unique<MTArrayT<MT>>(components, MT::Configuration
          (R, true, false, true, .5, 50));
      }
    }
    case 3: {
      typedef MonTableDivList<false,false> MT;
      return make_unique<MTArrayT<MT>>(components, MT::Configuration
        (R, true, false, true, .5, 500));
    }
    case 4: // fall through
    default: {
      if (allowRemovals) {
        typedef MonTableKDTree<false,false,1,true> MT;
        return make_unique<MTArrayT<MT>>(components, MT::Configuration
          (R, true, false, true, .5, 50));
      } else {
        typedef MonTableKDTree<false,false,1,false> MT;
        return make_unique<MTArrayT<MT>>(components, MT::Configuration
          (R, true, false, true, .5, 50));
      }
    }
  }
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
