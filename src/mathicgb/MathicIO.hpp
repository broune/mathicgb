#ifndef MATHICGB_MATHIC_IO_GUARD
#define MATHICGB_MATHIC_IO_GUARD

#include "Basis.hpp"
#include "Poly.hpp"
#include "Scanner.hpp"
#include "PolyRing.hpp"
#include "MonoProcessor.hpp"
#include <ostream>
#include <string>

/// Class for input and output in Mathic's format.
class MathicIO {
public:
  typedef PolyRing::Field BaseField;
  typedef BaseField::Element Coefficient;
  typedef BaseField::RawElement RawCoefficient;

  typedef PolyRing::Monoid Monoid;
  typedef Monoid::MonoRef MonoRef;
  typedef Monoid::Exponent Exponent;
  typedef Monoid::VarIndex VarIndex;
  typedef Monoid::ConstMonoRef ConstMonoRef;
  typedef Monoid::Order Order;
  typedef MonoProcessor<Monoid> Processor;

  typedef Order::Gradings Gradings;

/*
  /// reads ring, #gens, each generator in turn
  typedef std::tuple<
    std::unique_ptr<PolyRing>,
    std::unique_ptr<Basis>,
    std::unique_ptr<MonoProcessor<PolyRing::Monoid>>
  > BasisData;
  BasisData readBasis();



*/
  BaseField readBaseField(Scanner& in);
  void writeBaseField(const BaseField& field, std::ostream& out);

  std::pair<std::unique_ptr<PolyRing>, Processor> readRing(
    const bool withComponent,
    Scanner& in
  );
  void writeRing(
    const PolyRing& ring,
    const Processor& processor,
    const bool withComponent,
    std::ostream& out
  );

  Order readOrder(
    const VarIndex varCount,
    const bool withComponent,
    Scanner& in
  );

  void writeOrder(
    const Order& order,
    const bool withComponent,
    std::ostream& out
  );

  Basis readBasis(
    const PolyRing& ring,
    const bool readComponent,
    Scanner& in
  );

  void writeBasis(
    const Basis& basis,
    const bool writeComponent,
    std::ostream& out
  );

  Poly readPoly(
    const PolyRing& ring,
    const bool readComponent,
    Scanner& in
  );

  void writePoly(
    const Poly& poly,
    const bool writeComponent,
    std::ostream& out
  );

  void readTerm(
    const PolyRing& ring,
    const bool readComponent,
    Coefficient& coef,
    MonoRef mono,
    Scanner& in
  );

  void writeTerm(
    const PolyRing& ring,
    const bool writeComponent,
    const Coefficient coef,
    ConstMonoRef mono,
    std::ostream& out
  );

  /// Read a monomial with no coefficient. A 1 on its own is not
  /// considered a coefficient here - it is the monomial with all
  /// exponents zero and no coefficient.
  ///
  /// @todo: Eventually, pick up readComponent from Monoid::HasComponent.
  void readMonomial(
    const Monoid& monoid,
    const bool readComponent,
    MonoRef mono,
    Scanner& in
  );

  /// Print a monomial with no coefficient.
  void writeMonomial(
    const Monoid& monoid,
    const bool writeComponent,
    ConstMonoRef mono,
    std::ostream& out
  );

  /// Read the trailing indicator of the component of a module monomial.
  void readComponent(
    const Monoid& monoid,
    MonoRef mono,
    Scanner& in
  );

  void writeComponent(
    const Monoid& monoid,
    ConstMonoRef mono,
    std::ostream& out
  );
};

#endif
