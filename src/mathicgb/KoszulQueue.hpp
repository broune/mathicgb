#ifndef koszul_queue_guard
#define koszul_queue_guard

#include <memtailor.h>
#include <mathic.h>
#include "PolyRing.hpp"

class FreeModuleOrder;

class KoszulQueue
{
public:
  KoszulQueue(const FreeModuleOrder *order,
              memt::BufferPool &pool) :
    mQueue(Configuration(order, pool))
  {
  }

  const_monomial top() const {
    ASSERT(!empty());
    return mQueue.top();
  }

  monomial popRelease() {
    ASSERT(!empty());
    return mQueue.pop();
  }

  void pop() {
    ASSERT(!empty());
    mQueue.getConfiguration().getPool().free(popRelease().unsafeGetRepresentation());
  }

  void push(monomial sig) {
    //TODO: ASSERT(mPool.member(sig));
    mQueue.push(sig);
  }

  bool empty() const { return mQueue.empty(); }

  size_t size() const { return mQueue.size(); }

  size_t getMemoryUse() const {return mQueue.getMemoryUse();}

private:
  class Configuration
  {
  public:
    typedef monomial Entry;

    Configuration(const FreeModuleOrder *order,
                  memt::BufferPool &pool) :
      mOrder(order),
      mPool(pool)
    {
      ASSERT(mOrder != 0);
    }

    memt::BufferPool &getPool() { return mPool; }

    typedef int CompareResult; /* LT, EQ, GT */

    CompareResult compare(const Entry& a, const Entry& b) const;

    bool cmpLessThan(CompareResult r) const {return r == GT;}

    static const bool fastIndex = false;
    static const bool supportDeduplication = true;
    bool cmpEqual(CompareResult r) const {return r == EQ;}

    Entry deduplicate(Entry a, Entry b)
    {
      mPool.free(b.unsafeGetRepresentation());
      return a;
    }
  private:
    const FreeModuleOrder *mOrder;
    memt::BufferPool &mPool; // for signature monomials
  };

  mic::Heap<Configuration> mQueue;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
