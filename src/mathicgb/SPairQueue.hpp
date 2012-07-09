#ifndef _sPairQueue_h_
#define _sPairQueue_h_

#include <ostream>
class SPairGroup;

// Get concrete instances through a module order.
class SPairQueue {
public:
  typedef SPairGroup* Entry;

  virtual bool empty() const = 0;
  virtual void push(const Entry& entry) = 0;
  virtual Entry pop() = 0;
  virtual Entry top() const = 0;
  virtual void decreaseTop(const Entry& newValue) = 0;

  virtual void print(std::ostream& out) const = 0;
  virtual std::string getName() const = 0;
  virtual size_t getMemoryUse() const = 0;
  virtual size_t sumOfSizes() const = 0;
};

#endif
