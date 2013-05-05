#ifndef MATHICGB_MONO_ORDER_GUARD
#define MATHICGB_MONO_ORDER_GUARD

#include <vector>
#include <algorithm>

/// Class used to describe an monomial order and/or a module monomial
/// order. Use this class to construct a monoid. The monoid does the
/// actual comparisons.
///
/// You can extend the functionality for ordering offered here for
/// module mononials by pre- and post-processing of monomials. See
/// MonoProcessor.
template<class Weight>
class MonoOrder;

template<class W>
class MonoOrder {
public:
  typedef W Weight;
  typedef size_t VarIndex;
  typedef std::vector<Weight> Gradings;

  static const size_t ComponentAfterBaseOrder = static_cast<size_t>(-1);

  enum BaseOrder {
    /// Lexicographic order with x0 < x1 < ... < x_n.
    LexBaseOrder = 0,

    /// Reverse lexicographic order with x0 > x1 > ... > x_n.
    RevLexBaseOrder = 1
  };

  /// Same as MonoOrder(varCount, varOrder, gradings, componentBefore)
  /// where gradings has a single row of varCount 1's.
  MonoOrder(
    const VarIndex varCount,
    const BaseOrder baseOrder = RevLexBaseOrder,
    const size_t componentBefore = ComponentAfterBaseOrder,
    const bool componentsAscendingDesired = true,
    const bool schreyering = true
  ):
    mVarCount(varCount),
    mGradings
      (addComponentGrading(Gradings(varCount, 1), varCount, componentBefore)),
    mBaseOrder(baseOrder),
    mComponentGradingIndex(componentBefore),
    mComponentsAscendingDesired(componentsAscendingDesired),
    mSchreyering(schreyering)
  {}

  /// The specified base order is graded by the gradings matrix.
  ///
  /// The layout of the gradings matrix is row-major. For comparisons,
  /// the degree with respect to the first row is considered first,
  /// then the degree with respect to the second row and so on. The
  /// base order is used as a tie-breaker. The gradings vector can be
  /// empty. The order must be a monomial order - in particular, 1
  /// must be strictly less than all other monomials.
  ///
  /// For module monomials, the component is considered too. When the
  /// component is considered depens on componentBefore. If
  /// componentBefore == 0 then the component is considered before
  /// anything else. If componentBefore < gradingCount(), then the
  /// component is considered before the grading with index
  /// componentBefore and after the grading with index componentBefore
  /// - 1. If componentBefore == gradingCount(), then the component is
  /// considered after all gradings and before the base order. If
  /// componentBefore = ComponentAfterBaseOrder then the component is
  /// considered after everything else.
  ///
  /// The ordering represented by this object will be equivalent to
  /// the specified order, but it may be encoded differently. For
  /// example, if varCount == 0, then all orders are equivalent so all
  /// the other parameters are ignored and the encoding will be chosen
  /// to be the minimal-overhead one: ungraded lex with component
  /// considered last (do not depend on this particular choice - it
  /// may change). Therefore, you cannot count on any setting of this
  /// order to match what you passed in as a parameter.
  ///
  /// If the ordering you hav specified is not a monomial order, then
  /// the only guarantee in terms of encoding is that isMonomialOrder
  /// will return false.
  MonoOrder(
    const VarIndex varCount,
    Gradings&& gradings,
    const BaseOrder baseOrder = RevLexBaseOrder,
    const size_t componentBefore = ComponentAfterBaseOrder,
    const bool componentsAscendingDesired = true,
    const bool schreyering = true
  ):
    mVarCount(varCount),
    mGradings(std::move(gradings)),
    mBaseOrder(baseOrder),
    mComponentGradingIndex(componentBefore),
    mComponentsAscendingDesired(componentsAscendingDesired),
    mSchreyering(schreyering)
  {
#ifdef MATHCGB_DEBUG
    if (componentBefore != ComponentAfterBaseOrder) {
      for (VarIndex var = 0; var < varCount(); ++var) {
        MATHICGB_ASSERT(mGradings[var] == 0);
      }
    }
#endif
  }

  VarIndex varCount() const {return mVarCount;}

  VarIndex componentGradingIndex() const {return mComponentGradingIndex;}

  /// Returns the number of rows in the grading vector.
  size_t gradingCount() const {
    return varCount() == 0 ? 0 : mGradings.size() / varCount();
  }

  /// Returns the grading matrix in row-major layout.
  const Gradings& gradings() const {return mGradings;}

  /// Returns true if the grading matrix is a single row of 1's.
  bool isTotalDegree() const {
    if (varCount() == 0 || mGradings.size() != varCount())
      return false;
    for (VarIndex var = 0; var < varCount(); ++var)
      if (mGradings[var] != 1)
        return false;
    return true;
  }

  BaseOrder baseOrder() const {return mBaseOrder;}

  /// Returns true if the order is a monomial order. A monomial order
  /// is a total order on monomials where a>b => ac>bc for all
  /// monomials a,b,c and where the order is a well order. Only the
  /// well order property could currently fail. It is equivalent to
  /// stating that x>1 for all variables x.
  bool isMonomialOrder() const {
    for (VarIndex var = 0; var < varCount(); ++var) {
      // Check that x_var > 1.
      for (size_t grading = 0; ; ++grading) {
        if (grading == gradingCount) {
          // The column was entirely zero, so x_var > 1 if and only if the
          // base ordering is lex.
          if (baseOrder() != LexBaseOrder)
            return false;
          break;
        }
        const auto index = grading * gradingCount() + var;
        MATHICGB_ASSERT(index < mGradings.size());
        const auto weight = mGradings[index];
        if (weight != 0) {
          // We have found the first non-zero weight in this column,
          // so x_var > 1 if and only if this weight is positive.
          if (weight < 0)
            return false;
          break;
        }
      }
    }
    return true;
  }

  bool componentsAscendingDesired() const {return mComponentsAscendingDesired;}
  bool schreyering() const {return mSchreyering;}

private:
  static Gradings addComponentGrading(
    Gradings&& gradings,
    const VarIndex varCount,
    const VarIndex componentBefore
  ) {
    if (componentBefore == ComponentAfterBaseOrder)
      return std::move(gradings);
    MATHICGB_ASSERT(componentBefore <= varCount);
    gradings.resize(gradings.size() + varCount);
    const auto newRow = gradings.begin() + varCount * componentBefore;
    std::copy_n(newRow, varCount, newRow + varCount);
    std::fill_n(newRow, varCount, static_cast<Weight>(0));
    return std::move(gradings);
  }

  bool debugAssertValid() {
#ifdef MATHICGB_DEBUG
    MATHICGB_ASSERT(mGradings.size() == gradingCount() * varCount());
    MATHICGB_ASSERT(
      mComponentGradingIndex == ComponentAfterBaseOrder ||
      mComponentGradingIndex < gradingCount()
    );
    MATHICGB_ASSERT(
      mBaseOrder == LexBaseOrder ||
      mBaseOrder == RevLexBaseOrder
    );
    if (varCount() == 0) {
      MATHICGB_ASSERT(mGradings.empty());
      MATHICGB_ASSERT(baseOrder() == RevLexBaseOrder());
      MATHICGB_ASSERT(mComponentGradingIndex == ComponentAfterBaseOrder);
    }
#endif
    return true;
  }

  const VarIndex mVarCount;

  const Gradings mGradings;
  const BaseOrder mBaseOrder;
  const size_t mComponentGradingIndex;
  const bool mSchreyering;
  const bool mComponentsAscendingDesired;
};

#endif
