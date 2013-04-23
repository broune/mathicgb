#ifndef MATHICGB_MONO_ORDER_GUARD
#define MATHICGB_MONO_ORDER_GUARD

#include <vector>

/// Class used to describe an monomial order or a module monomial
/// order. Use this class to construct a monoid. The monoid does the
/// actual comparisons.
///
/// You can extend the functionality for ordering offered here for
/// module mononials by pre- and post-processing. See MonoProcessor.
template<class Weight>
class MonoOrder;

template<class W>
class MonoOrder {
public:
  typedef W Weight;
  typedef size_t VarIndex;

  MonoOrder(const VarIndex varCount):
    mVarCount(varCount),
    mGradings(false),
    mTotalDegreeGrading(true),
    mBaseOrder(RevLexBaseOrder),
    mComponentBefore(ComponentLast),
    mUseSchreyerOrder(true),
    mComponentsAscendingDesired(true)
  {}

  const VarIndex varCount() const {return mVarCount;}

  /// Set the matrix that defines the order. The entry layout is row
  /// major. The first row is considered first, then the second row
  /// and so on.
  void setGradings(std::vector<Weight>&& gradings) {
    MATHICGB_ASSERT(varCount() > 0 || gradings.empty());
    MATHICGB_ASSERT(varCount() == 0 || gradings.size() % varCount() == 0);
    mGradings = std::move(gradings);
    mTotalDegreeGrading = false;
  }

  /// As calling setGradings with a vector of varCount() 1's.
  void setTotalDegreeGrading() const {
    mGradings.clear();
    mTotalDegreeGrading = true;
  }

  std::vector<Weight>& gradings() {
    if (mTotalDegreeGrading && mGradings.empty())
      mGradings.resize(varCount(), 1);
    return mGradings;
  }

  bool gradingIsTotalDegreeRevLex() const {
    return baseOrder() == RevLexBaseOrder && mTotalDegreeGrading;
  }

  /// Returns the number of rows in the grading vector.
  size_t gradingCount() const {return mGradings.size() / varCount();}

  enum BaseOrder {
    LexBaseOrder = 0,
    RevLexBaseOrder = 1
  };
  /// Set the order to use as a tie-braker for the grading matrix.
  void setBaseOrder(BaseOrder baseOrder) {mBaseOrder = baseOrder;}
  BaseOrder baseOrder() const {return mBaseOrder;}

  static const size_t ComponentLast = static_cast<size_t>(-1);

  /// The component is considered before row (grading) and after row
  /// (grading - 1) in the order matrix. If grading == 0 then the
  /// component is considered before anything else. If grading ==
  /// gradingCount() then grading is considered after the grading
  /// matrix but before the base order. If grading is ComponentLast,
  /// then the component is used as the final tie breaker after both
  /// the grading matrix and the base order. All values of this setting
  /// except for ComponentLast imply overhead equivalent to one extra row
  /// in the grading matrix even for monoids without components.
  /// ComponentLast has no overhead.
  void setComponentBefore(const size_t grading) const {
    MATHICGB_ASSERT(grading == ComponentLast || grading <= gradingCount);
    mComponentBefore = grading;
  }

private:
  const VarIndex mVarCount;
  std::vector<Weight> mGradings;
  bool mTotalDegreeGrading;
  BaseOrder mBaseOrder;
  size_t mComponentBefore;
  bool mUseSchreyerOrder;
  bool mComponentsAscendingDesired;
};

#endif
