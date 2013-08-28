// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#ifndef MATHICGB_M_T_ARRAY_GUARD
#define MATHICGB_M_T_ARRAY_GUARD

#include "PolyRing.hpp"
#include <vector>
#include <iostream>
#include <algorithm>

MATHICGB_NAMESPACE_BEGIN

class MonomialTableArray
{
public:
  virtual ~MonomialTableArray() {};

  // returns true if the monomial actually needs to be inserted.
  // If the monomial is inserted, the caller agrees to keep that monomial
  // active until it is removed from the table.
  virtual bool insert(const_monomial m) = 0;

  virtual bool member(const_monomial m) = 0;

  virtual void display(std::ostream& o) const = 0;

  virtual void getMonomials(std::vector<const_monomial>& monomials) const = 0;

  virtual size_t n_elems() const = 0;

  virtual std::string name() const = 0;

  virtual size_t getMemoryUse() const = 0;

  // Choosing one
  static int displayMTTypes(std::ostream &o); // returns n s.t. 0..n-1 are valid types
  static std::unique_ptr<MonomialTableArray>
    make(const PolyRing *R, int typ, size_t components, bool allowRemovals);
};

MATHICGB_NAMESPACE_END
#endif
