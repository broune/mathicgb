// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_ZIP_GUARD
#define MATHICGB_ZIP_GUARD

#include "Range.hpp"

MATHICGB_NAMESPACE_BEGIN

/// Zip ties two iterators together into a single iterator. The point
/// is to enable the zip() function defined further down in this header.
template<class Iterator1, class Iterator2>
class Zip {
public:
  typedef decltype(*std::declval<Iterator1>()) Value1;
  typedef decltype(*std::declval<Iterator2>()) Value2;
  typedef std::pair<Value1, Value2> value_type;

  Zip() {}
  Zip(Iterator1 it1, Iterator2 it2): mIt1(it1), mIt2(it2) {}

  Zip& operator++() {
    ++mIt1;
    ++mIt2;
    return *this;
  }

  value_type operator*() const {
    // We cannot use std::make_pair here. If .first or .second is a reference,
    // then std::make_pair will remove the reference and then the conversion
    // to value_type will add the reference back in - but now the reference
    // is not what was originally referenced, instead it is a reference to
    // the temporary object returned by std::make_pair, which then promptly
    // goes out of scope, leading to undefined behavior.
    return value_type(*mIt1, *mIt2);
  }

  bool operator!=(const Zip& it) const {
    MATHICGB_ASSERT((mIt1 != it.mIt1) == (mIt2 != it.mIt2));
    return mIt1 != it.mIt1;
  }

  bool operator==(const Zip& it) const {
    MATHICGB_ASSERT((mIt1 == it.mIt1) == (mIt2 == it.mIt2));
    return mIt1 == it.mIt1;
  }

private:
  Iterator1 mIt1;
  Iterator2 mIt2;
};

/// Creates a Zip iterator out of it1 and it2. This is a convenience function
/// that avoids the need to specify the iterator types explicitly.
template<class Iterator1, class Iterator2>
auto makeZip(
  Iterator1&& it1,
  Iterator2&& it2
) -> Zip<Iterator1, Iterator2> {
  return Zip<Iterator1, Iterator2>(
    std::forward<Iterator1>(it1),
    std::forward<Iterator2>(it2)
  );
}

/// Zips two ranges into a single zipped range. Example:
///
///   std::vector<string> a = {"hello", "world"};
///   std::vector<int> b = {4, 2};
///   for (const auto& p : zip(a, b))
///     std::cout << p.first << p.second << ' ';
///
/// The output will be "hello4 world2 ".
template<class Range1, class Range2>
auto zip(Range1&& range1, Range2&& range2) ->
  decltype(
    range(
      makeZip(std::begin(range1), std::begin(range2)),
      makeZip(std::end(range1), std::end(range2))
    )
  )
{
  return range(
    makeZip(std::begin(range1), std::begin(range2)),
    makeZip(std::end(range1), std::end(range2))
  );
}

MATHICGB_NAMESPACE_END

#endif
