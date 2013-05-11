#ifndef koszul_queue_guard
#define koszul_queue_guard

#include "PolyRing.hpp"
#include "NonCopyable.hpp"
#include <mathic.h>

class KoszulQueue : public NonCopyable<KoszulQueue> {
public:
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::Mono Mono;
  typedef Monoid::ConstMonoRef ConstMonoRef;

  KoszulQueue(const Monoid& monoid): mQueue(Configuration(monoid)) {}
  KoszulQueue(KoszulQueue&& kq): mQueue(std::move(kq.mQueue)) {}

  ConstMonoRef top() const {
    MATHICGB_ASSERT(!empty());
    return mQueue.top();
  }

  void pop() {
    MATHICGB_ASSERT(!empty());
    monoid().freeRaw(mQueue.pop());
  }

  void push(ConstMonoRef sig) {
    auto m = monoid().alloc();
    monoid().copy(sig, m);
    mQueue.push(*m.release());
  }
  bool empty() const {return mQueue.empty();}
  size_t size() const {return mQueue.size();}

  size_t getMemoryUse() const {return mQueue.getMemoryUse();}

  const Monoid& monoid() const {return mQueue.getConfiguration().monoid();}

private:
  class Configuration {
  public:
    typedef Monoid::MonoRef Entry;

    Configuration(const Monoid& monoid): mMonoid(monoid) {}

    typedef Monoid::CompareResult CompareResult;
    CompareResult compare(const Entry& a, const Entry& b) const {
      return ring().monoid().compare(a, b);
    }
    bool cmpLessThan(CompareResult r) const {return r == Monoid::GreaterThan;}

    static const bool fastIndex = false;
    static const bool supportDeduplication = true;
    bool cmpEqual(CompareResult r) const {return r == Monoid::EqualTo;}

    Entry deduplicate(Entry a, Entry b)
    {
      monoid().freeRaw(b);
      return a;
    }

    const Monoid& monoid() const {return mMonoid;}

  private:
    const Monoid& mMonoid;
  };

  mic::Heap<Configuration> mQueue;
};

#endif

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
