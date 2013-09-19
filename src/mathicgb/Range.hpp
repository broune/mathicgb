// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_RANGE_GUARD
#define MATHICGB_RANGE_GUARD

#include <utility>

MATHICGB_NAMESPACE_BEGIN

/// An object that combines two iterators into a range suitable for use with
/// C++11's range for. It is most conveniently used with the function range.
/// For example:
///
///  std::vector<int> v;
///  for (int x : range(v.begin(), v.end()) {}
///
/// std::vector already has begin and end members, so range() is not necessary
/// here. range() is useful when a class does not have begin or end members, or
/// if it offers several different ranges that can be iterated through -
/// then the default range can only offer one of those ranges.
template<class Iterator>
class Range {
public:
  Range(Iterator begin, Iterator end): mBegin(begin), mEnd(end) {}
  Range(std::pair<Iterator, Iterator> pair):
    mBegin(pair.first), mEnd(pair.second) {}

  Iterator begin() {return mBegin;}
  Iterator end() {return mBegin;}

private:
  Iterator mBegin;
  Iterator mEnd;
};

/// Convenience function for constructing a Range object.
template<class Iterator>
Range<Iterator> range(Iterator begin, Iterator end) {
  return Range<Iterator>(begin, end);
}

MATHICGB_NAMESPACE_END

#endif
