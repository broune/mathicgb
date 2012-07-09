#include "stdinc.h"
#include "BitTriangle.hpp"

size_t BitTriangle::getMemoryUse() const
{
  size_t sum = mColumns.capacity() * sizeof(mColumns.front());
  const size_t stop = mColumns.size();
  for (size_t i = 0; i != stop; ++i) {
    size_t const capacity = mColumns[i].capacity();
    sum += (capacity / 8) + (capacity % 8 != 0); // 8 bits per byte rounded up
  }
  return sum;
}
