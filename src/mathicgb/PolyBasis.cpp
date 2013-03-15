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

  if (mUseBuchbergerLcmHitCache)
    mBuchbergerLcmHitCache.push_back(0);
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

bool PolyBasis::buchbergerLcmCriterion(size_t a, size_t b) const {
  MATHICGB_ASSERT(a != b);
  MATHICGB_ASSERT(!retired(a));
  MATHICGB_ASSERT(!retired(b));
  // We don't need to set the weights on the lcm since we will only use it
  // for testing divisibility, not for determining order.
  monomial lcmAB = mRing.allocMonomial();
  mRing.monomialLeastCommonMultipleNoWeights
    (leadMonomial(a), leadMonomial(b), lcmAB);
  bool value = buchbergerLcmCriterion(a, b, lcmAB);
  mRing.freeMonomial(lcmAB);
  return value;
}

bool PolyBasis::buchbergerLcmCriterion
  (size_t a, size_t b, const_monomial lcmAB) const
{
  MATHICGB_ASSERT(a != b);
  MATHICGB_ASSERT(!retired(a));
  MATHICGB_ASSERT(!retired(b));
  MATHICGB_ASSERT(mRing.monomialIsLeastCommonMultipleNoWeights
    (leadMonomial(a), leadMonomial(b), lcmAB));

  class Criterion : public DivisorLookup::EntryOutput {
  public:
    Criterion(size_t a, size_t b, const_monomial lcmAB, const PolyBasis& basis):
      mA(a), mB(b),
      mLcmAB(lcmAB),
      mRing(basis.ring()),
      mBasis(basis),
      mHit(static_cast<size_t>(-1)),
      mAlmostApplies(false) {}

    virtual bool proceed(size_t index) {
      MATHICGB_ASSERT(index < mBasis.size());
      MATHICGB_ASSERT(!applies()); // should have stopped search in this case
      MATHICGB_ASSERT(
        mRing.monomialIsDivisibleBy(mLcmAB, mBasis.leadMonomial(index)));
      if (!mBasis.leadMinimal(index))
        return true;
      if (index == mA || index == mB)
        return true;
      mAlmostApplies = true;

      if (index < mA || index < mB) {
        // check lcm(a,index) != lcm(a,b) <=>
        // exists i such that max(a[i], c[i]) != max(a[i],b[i]) <=>
        // exists i such that b[i] > a[i] && b[i] > c[i]
        const_monomial leadA = mBasis.leadMonomial(mA);
        const_monomial leadB = mBasis.leadMonomial(mB);
        const_monomial leadC = mBasis.leadMonomial(index);
        if (//index > mA &&
          !mRing.monomialHasStrictlyLargerExponent(leadB, leadC, leadA))
          return true;

        // check lcm(b,index) != lcm(a,b)
        if (//index > mB &&
          !mRing.monomialHasStrictlyLargerExponent(leadA, leadC, leadB))
          return true;
      }

      mHit = index;
      return false; // stop search
    }

    const_monomial lcmAB() const {return mLcmAB;}
    bool almostApplies() const {return mAlmostApplies;}
    bool applies() const {return mHit != static_cast<size_t>(-1);}
    size_t hit() const {return mHit;}

  private:
    size_t mA;
    size_t mB;
    const_monomial mLcmAB;
    const PolyRing& mRing;
    const PolyBasis& mBasis;
    size_t mHit; // the divisor that made the criterion apply
    bool mAlmostApplies; // applies ignoring lcm(a,b)=lcm(a,c) complication
  };

  ++mStats.buchbergerLcmQueries;
  bool applies = false;
  bool almostApplies = false;
  {
    Criterion criterion(a, b, lcmAB, *this);
    if (mUseBuchbergerLcmHitCache) {
      // Check cacheB first since when I tried this there was a higher hit rate
      // for cacheB than cacheA. Might not be a persistent phenomenon, but
      // there's no downside to trying out cacheB first so I'm going for that.
      //
      // I update the cache if the second check is a hit but not if the first
      // check is a hit. In the one test I did, the worst hit rate was from
      // updating the cache every time, the second best hit rate was from
      // not updating the cache (from cache hits) and the best hit rate was
      // from doing this.
      //
      // The idea is that when the first cache check is a hit,
      // the second cache member might have been a hit too, and updating it
      // might replace a high hit rate element with a low hit rate element,
      // which would be bad. When the second cache check is a hit, we know
      // that the first one wasn't (or we would have taken an early exit),
      // so we have reason to suspect that the first cache element is not
      // a high hit rate element. So it should be better to replace it.
      // That idea seems to be right since it worked better in the one
      // test I did.
      size_t cacheB = mBuchbergerLcmHitCache[b];
      if (!applies && !retired(cacheB) &&
        ring().monomialIsDivisibleBy
          (criterion.lcmAB(), leadMonomial(cacheB)))
        applies = !criterion.Criterion::proceed(cacheB);

      size_t cacheA = mBuchbergerLcmHitCache[a];
      if (!applies && !retired(cacheA) &&
        ring().monomialIsDivisibleBy
          (criterion.lcmAB(), leadMonomial(cacheA))) {
        applies = !criterion.Criterion::proceed(cacheA);
        if (applies)
          mBuchbergerLcmHitCache[b] = cacheA;
      }
    }
    if (applies)
      ++mStats.buchbergerLcmCacheHits;
    else {
      MATHICGB_ASSERT(!criterion.applies());
      mDivisorLookup->divisors(criterion.lcmAB(), criterion);
      applies = criterion.applies();

      if (mUseBuchbergerLcmHitCache && applies) {
        MATHICGB_ASSERT(criterion.hit() < size());
        mBuchbergerLcmHitCache[a] = criterion.hit();
        mBuchbergerLcmHitCache[b] = criterion.hit();
      }
    }
    if (!applies)
      almostApplies = criterion.almostApplies();
  }

  if (applies)
    ++mStats.buchbergerLcmHits;
  else if (almostApplies)
    ++mStats.buchbergerLcmNearHits;

  MATHICGB_ASSERT(applies == buchbergerLcmCriterionSlow(a, b));
  return applies;
}

bool PolyBasis::buchbergerLcmCriterionSlow(size_t a, size_t b) const {
  MATHICGB_ASSERT(a != b);
  MATHICGB_ASSERT(!retired(a));
  MATHICGB_ASSERT(!retired(b));

  monomial lcmAB = ring().allocMonomial();
  monomial lcm = ring().allocMonomial();
  ring().monomialLeastCommonMultiple
    (leadMonomial(a), leadMonomial(b), lcmAB);
  size_t stop = size();
  size_t i = 0;
  for (; i < stop; ++i) {
    if (retired(i) || !leadMinimal(i))
      continue;
    if (!ring().monomialIsDivisibleBy(lcmAB, leadMonomial(i)))
      continue;
    if (i == a || i == b)
      continue;
    if (i < a || i < b) {
      ring().monomialLeastCommonMultiple
        (leadMonomial(a), leadMonomial(i), lcm);
      if (ring().monomialEQ(lcmAB, lcm))
        continue;

      ring().monomialLeastCommonMultiple
        (leadMonomial(b), leadMonomial(i), lcm);
      if (ring().monomialEQ(lcmAB, lcm))
        continue;
    }
    break;
  }
  ring().freeMonomial(lcmAB);
  ring().freeMonomial(lcm);
  return i != stop;
}

PolyBasis::Entry::Entry():
  poly(0),
  leadMinimal(0),
  retired(false),
  usedAsStartCount(0),
  usedAsReducerCount(0),
  possibleReducerCount(0),
  nonSignatureReducerCount(0) {}
