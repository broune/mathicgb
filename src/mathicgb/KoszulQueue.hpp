#ifndef koszul_queue_guard
#define koszul_queue_guard

#include "PolyRing.hpp"
#include "NonCopyable.hpp"
#include <mathic.h>

class KoszulQueue : public NonCopyable<KoszulQueue> {
public:
  KoszulQueue(const PolyRing& ring): mQueue(Configuration(ring)) {}
  KoszulQueue(KoszulQueue&& kq): mQueue(std::move(kq.mQueue)) {}

  const_monomial top() const {
    MATHICGB_ASSERT(!empty());
    return mQueue.top();
  }

  monomial popRelease() {
    MATHICGB_ASSERT(!empty());
    return mQueue.pop();
  }

  void pop() {
    MATHICGB_ASSERT(!empty());
    auto toFree = popRelease().unsafeGetRepresentation();
    mQueue.getConfiguration().ring().freeMonomial(toFree);
  }

  void push(monomial sig) {mQueue.push(sig);}
  bool empty() const {return mQueue.empty();}
  size_t size() const {return mQueue.size();}

  size_t getMemoryUse() const {return mQueue.getMemoryUse();}

private:
  KoszulQueue(const KoszulQueue&); // unavailable
  

  class Configuration
  {
  public:
    typedef monomial Entry;

    Configuration(const PolyRing& ring): mRing(ring) {}

    const PolyRing& ring() const {return mRing;}

    typedef int CompareResult; /* LT, EQ, GT */

    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monoid().compare(a, b);
    }

    bool cmpLessThan(CompareResult r) const {return r == GT;}

    static const bool fastIndex = false;
    static const bool supportDeduplication = true;
    bool cmpEqual(CompareResult r) const {return r == EQ;}

    Entry deduplicate(Entry a, Entry b)
    {
      ring().freeMonomial(b.unsafeGetRepresentation());
      return a;
    }
  private:
    const PolyRing& mRing;
  };

  mic::Heap<Configuration> mQueue;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
