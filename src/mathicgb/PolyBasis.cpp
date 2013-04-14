#include "stdinc.h"
#include "Ideal.hpp"
#include "PolyBasis.hpp"

#include "DivisorLookup.hpp"

PolyBasis::PolyBasis(
  const PolyRing& ring,
  FreeModuleOrder& order,
  std::unique_ptr<DivisorLookup> divisorLookup
):
  mRing(ring),
  mOrder(order),
  mDivisorLookup(std::move(divisorLookup))
{
  MATHICGB_ASSERT(mDivisorLookup.get() != 0);
  mDivisorLookup->setBasis(*this);
}

PolyBasis::~PolyBasis() {
  EntryIter const stop = mEntries.end();
  for (EntryIter it = mEntries.begin(); it != stop; ++it) {
    if (it->retired)
      continue;
    MATHICGB_ASSERT(it->poly != 0);
    delete it->poly;
  }
}

std::unique_ptr<Ideal> PolyBasis::initialIdeal() const {
  std::unique_ptr<Ideal> ideal(new Ideal(mRing));
  size_t const idealSize = size();
  for (size_t gen = 0; gen != idealSize; ++gen) {
    if (!retired(gen) && leadMinimal(gen)) {
      std::unique_ptr<Poly> p(new Poly(mRing));
      p->appendTerm(1, leadMonomial(gen));
      ideal->insert(std::move(p));
    }
  }
  ideal->sort(mOrder);
  return ideal;
}

void PolyBasis::insert(std::unique_ptr<Poly> poly) {
  MATHICGB_ASSERT(poly.get() != 0);
  MATHICGB_ASSERT(!poly->isZero());
  poly->makeMonic();
  const size_t index = size();
  EntryIter const stop = mEntries.end();
  const_monomial const lead = poly->getLeadMonomial();
#ifdef DEBUG
  // lead monomials must be unique among basis elements
  for (EntryIter it = mEntries.begin(); it != stop; ++it) {
    if (it->retired)
      continue;
    MATHICGB_ASSERT(!ring().monomialEQ(lead, it->poly->getLeadMonomial()));
  }
#endif

  // Update information about minimal lead terms.
  const bool leadMinimal = (divisor(lead) == static_cast<size_t>(-1));
  if (leadMinimal) {
    class MultipleOutput : public DivisorLookup::EntryOutput {
    public:
      MultipleOutput(EntryCont& entries): mEntries(entries) {}
      virtual bool proceed(size_t index) {
        MATHICGB_ASSERT(index < mEntries.size());
        mEntries[index].leadMinimal = false;
        return true;
      }
    private:
      EntryCont& mEntries;
    };
    MultipleOutput out(mEntries);
    divisorLookup().multiples(lead, out);
  }

  mDivisorLookup->insert(lead, index);

  mEntries.push_back(Entry());
  Entry& entry = mEntries.back();
  entry.poly = poly.release();
  entry.leadMinimal = leadMinimal;

  MATHICGB_ASSERT(mEntries.back().poly != 0);
}

std::unique_ptr<Poly> PolyBasis::retire(size_t index) {
  MATHICGB_ASSERT(index < size());
  MATHICGB_ASSERT(!retired(index));
  mDivisorLookup->remove(leadMonomial(index));
  std::unique_ptr<Poly> poly(mEntries[index].poly);
  mEntries[index].poly = 0;
  mEntries[index].retired = true;
  return poly;
}

std::unique_ptr<Ideal> PolyBasis::toIdealAndRetireAll() {
  auto ideal = make_unique<Ideal>(ring());
  for (size_t i = 0; i < size(); ++i)
    if (!retired(i))
      ideal->insert(retire(i));
  return ideal;
}

size_t PolyBasis::divisor(const_monomial mon) const {
  size_t index = divisorLookup().divisor(mon);
  MATHICGB_ASSERT((index == static_cast<size_t>(-1)) ==
    (divisorSlow(mon) == static_cast<size_t>(-1)));
  MATHICGB_ASSERT(index == static_cast<size_t>(-1) ||
    ring().monomialIsDivisibleBy(mon, leadMonomial(index)));
  return index;
}

size_t PolyBasis::classicReducer(const_monomial mon) const {
  return divisorLookup().classicReducer(mon);
  size_t index = divisorLookup().classicReducer(mon);
  MATHICGB_ASSERT((index == static_cast<size_t>(-1)) ==
    (divisorSlow(mon) == static_cast<size_t>(-1)));
  MATHICGB_ASSERT(index == static_cast<size_t>(-1) ||
    ring().monomialIsDivisibleBy(mon, leadMonomial(index)));
  return index;
}

size_t PolyBasis::divisorSlow(const_monomial mon) const {
  const size_t stop = size();
  for (size_t i = 0; i != stop; ++i)
    if (!retired(i) && ring().monomialIsDivisibleBy(mon, leadMonomial(i)))
      return i;
  return static_cast<size_t>(-1);
}

bool PolyBasis::leadMinimalSlow(size_t index) const {
  MATHICGB_ASSERT(index < size());
  MATHICGB_ASSERT(!retired(index));
  const_monomial const lead = leadMonomial(index);
  EntryCIter const skip = mEntries.begin() + index;
  EntryCIter const stop = mEntries.end();
  for (EntryCIter it = mEntries.begin(); it != stop; ++it) {
    if (it->retired)
      continue;
    const_monomial const itLead = it->poly->getLeadMonomial();
    if (ring().monomialIsDivisibleBy(lead, itLead) && it != skip)
      return false;
  }
  return true;
}

size_t PolyBasis::minimalLeadCount() const {
  // todo: use iterators
  size_t minCount = 0;
  const size_t stop = size();
  for (size_t i = 0; i != stop; ++i)
    if (!retired(i) && leadMinimal(i))
      ++minCount;
  return minCount;
}

size_t PolyBasis::maxIndexMinimalLead() const {
  // todo: use iterators
  size_t i = size() - 1;
  for (; i != static_cast<size_t>(-1); --i)
    if (!retired(i) && leadMinimal(i))
      break;
  return i;
}

size_t PolyBasis::monomialCount() const {
  size_t sum = 0;
  EntryCIter const stop = mEntries.end();
  for (EntryCIter it = mEntries.begin(); it != stop; ++it)
    if (!it->retired)
      sum += it->poly->nTerms();
  return sum;
}

size_t PolyBasis::getMemoryUse() const {
  size_t sum = mEntries.capacity() * sizeof(mEntries.front());
  EntryCIter const stop = mEntries.end();
  for (EntryCIter it = mEntries.begin(); it != stop; ++it)
    if (!it->retired)
      sum += it->poly->getMemoryUse();
  return sum;
}

PolyBasis::Entry::Entry():
  poly(0),
  leadMinimal(0),
  retired(false),
  usedAsStartCount(0),
  usedAsReducerCount(0),
  possibleReducerCount(0),
  nonSignatureReducerCount(0) {}
