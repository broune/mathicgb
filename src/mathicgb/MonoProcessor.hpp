#ifndef MATHICGB_MONO_PROCESSOR_GUARD
#define MATHICGB_MONO_PROCESSOR_GUARD

#include "MonoMonoid.hpp"

/// Does pre- and post-processing of monomials to implement monomial
/// orders not directly supported by the monoid. This is so far only
/// relevant for module monomials.
///
/// todo: distinguish monomials from module monomials using two
/// different monoids.
template<class Monoid>
class MonoProcessor;

template<class M>
class MonoProcessor {
public:
  typedef M Monoid;
  typedef typename Monoid::VarIndex VarIndex;
  typedef typename Monoid::MonoVector MonoVector;
  typedef typename Monoid::MonoRef MonoRef;
  typedef typename Monoid::ConstMonoRef ConstMonoRef;
  typedef typename Monoid::ConstMonoPtr ConstMonoPtr;

  static_assert(Monoid::HasComponent, "");

  MonoProcessor(
    bool componentsAscendingDesired, // todo: replace with MonoOrder
    VarIndex componentCount, // todo: replace with MonoOrder
    const Monoid& monoid,
    MonoVector&& moduleAdjustments
  ):
    mMonoid(monoid),
    mComponentsAscendingDesired(componentsAscendingDesired),
    mComponentCount(componentCount),
    mModuleAdjustmentsMemory(std::move(moduleAdjustments))
  {
    MATHICGB_ASSERT(mModuleAdjustmentsMemory.monoid() == this->monoid());
    MATHICGB_ASSERT(mModuleAdjustmentsMemory.empty() ||
      mModuleAdjustmentsMemory.size() == this->componentCount());
    // todo: add a component count to the monoid and get it from there
    // instead of storing it here.

    for (
      auto it = mModuleAdjustmentsMemory.begin();
      it != mModuleAdjustmentsMemory.end();
      ++it
    ) {
      // in the absence of a separate monoid for (non-module) monomials,
      // at least we can check that the component is zero.
      MATHICGB_ASSERT(this->monoid().component(*it) == 0);

      // todo: there should be a better way of indexing into a
      // MonoVector.
      mModuleAdjustments.emplace_back((*it).ptr());
    }
  }
    

  void preprocess(MonoRef mono) const {
    if (hasModuleAdjustments())
      monoid().multiplyInPlace(moduleAdjustment(mono), mono);
    if (needToReverseComponents())
      reverseComponent(mono);
  }

  void postprocess(MonoRef mono) const {
    if (needToReverseComponents())
      reverseComponent(mono);
    if (hasModuleAdjustments()) {
      MATHICGB_ASSERT(monoid().divides(moduleAdjustment(mono), mono));
      monoid().divideInPlace(moduleAdjustment(mono), mono);
    }
  }

  bool processingNeeded() const {
    return needToReverseComponents() || hasModuleAdjustments();
  }

  bool needToReverseComponents() const {
    return Monoid::HasComponent &&
      mComponentsAscendingDesired != mMonoid.componentsAscending();
  }

  bool hasModuleAdjustments() const {
    return !mModuleAdjustments.empty();
  }

  VarIndex componentCount() const {return mComponentCount;}
  const Monoid& monoid() const {return mMonoid;}

private:
  void operator==(const MonoProcessor&) const; // not available

  void reverseComponent(MonoRef mono) const {
    const auto component = monoid().component(mono);
    const auto newComponent = mComponentCount - 1 - component;
    monoid().setComponent(newComponent, mono);
  }

  ConstMonoRef moduleAdjustment(ConstMonoRef mono) const {
    MATHICGB_ASSERT(hasModuleAdjustments());
    const auto component = monoid().component(mono);
    MATHICGB_ASSERT(component < componentCount());
    MATHICGB_ASSERT(mModuleAdjustments.size() == componentCount());
    return *mModuleAdjustments[component];
  }

  const Monoid& mMonoid;
  const bool mComponentsAscendingDesired;
  const VarIndex mComponentCount;
  MonoVector mModuleAdjustmentsMemory;
  std::vector<ConstMonoPtr> mModuleAdjustments;
};

#endif
