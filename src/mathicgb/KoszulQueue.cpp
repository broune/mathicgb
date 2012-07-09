#include "stdinc.h"
#include "KoszulQueue.hpp"
#include "FreeModuleOrder.hpp"

KoszulQueue::Configuration::CompareResult KoszulQueue::Configuration::compare(const Entry& a, const Entry& b) const {
  return mOrder->signatureCompare(a,b);
}

// Local Variables:
// compile-command: "make -C .. "
// indent-tabs-mode: nil
// End:
