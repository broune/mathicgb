#ifndef MATHICGB_MONO_PROCESSOR_GUARD
#define MATHICGB_MONO_PROCESSOR_GUARD

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

  MonoProcessor(const Monoid& monoid):
    mComponentsAscendingDesired(true),
    mComponentCount(0),
    mModuleAdjustmentsMemory(monoid)
  {}

  void setModuleAdjustments(MonoVector&& moduleAdjustments) {
    MATHICGB_ASSERT(moduleAdjustments.monoid() == monoid());
    MATHICGB_ASSERT(mModuleAdjustmentsMemory.empty() ||
      mModuleAdjustmentsMemory.size() == componentCount());
    mModuleAdjustmentsMemory = std::move(moduleAdjustments);

    mModuleAdjustments.clear();
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
      componentsAscendingDesired() != monoid().componentsAscending();
  }

  void setComponentsAscendingDesired(bool value) {
    mComponentsAscendingDesired = value;
  }
  bool componentsAscendingDesired() const {return mComponentsAscendingDesired;}

  bool hasModuleAdjustments() const {
    return !mModuleAdjustments.empty();
  }

  void setComponentCount(VarIndex count) {mComponentCount = count;}
  VarIndex componentCount() const {return mComponentCount;}
  const Monoid& monoid() const {return mModuleAdjustmentsMemory.monoid();}

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

  bool mComponentsAscendingDesired;
  VarIndex mComponentCount;
  MonoVector mModuleAdjustmentsMemory;
  std::vector<ConstMonoPtr> mModuleAdjustments;
};

#endif
