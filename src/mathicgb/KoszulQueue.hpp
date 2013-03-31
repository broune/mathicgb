#ifndef koszul_queue_guard
#define koszul_queue_guard

#include <memtailor.h>
#include <mathic.h>
#include "PolyRing.hpp"

class FreeModuleOrder;

class KoszulQueue
{
public:
  KoszulQueue(const FreeModuleOrder *order, const PolyRing& ring):
    mQueue(Configuration(order, ring))
  {}

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
  class Configuration
  {
  public:
    typedef monomial Entry;

    Configuration(const FreeModuleOrder *order, const PolyRing& ring):
      mOrder(order), mRing(ring)
    {
      MATHICGB_ASSERT(mOrder != 0);
    }

    const PolyRing& ring() const {return mRing;}

    typedef int CompareResult; /* LT, EQ, GT */

    CompareResult compare(const Entry& a, const Entry& b) const;

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
    const FreeModuleOrder *mOrder;
    const PolyRing& mRing;
  };

  mic::Heap<Configuration> mQueue;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
