#ifndef MATHICGB_MATHIC_IO_GUARD
#define MATHICGB_MATHIC_IO_GUARD

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

  std::pair<PolyRing, Processor> readRing(
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

  void readTerm(
    const PolyRing& ring,
    const bool readComponent,
    Coefficient& coef,
    MonoRef mono,
    Scanner& in
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

  /// Read the trailing indicator of the component of a module monomial.
  void readComponent(
    const Monoid& monoid,
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
};

auto MathicIO::readBaseField(Scanner& in) -> BaseField {
  return BaseField(in.readInteger<RawCoefficient>());
}

void MathicIO::writeBaseField(const BaseField& field, std::ostream& out) {
  out << field.charac();
}

auto MathicIO::readRing(
  const bool withComponent,
  Scanner& in
) -> std::pair<PolyRing, Processor> {
  auto baseField = readBaseField(in);
  const auto varCount = in.readInteger<VarIndex>();
  auto order = readOrder(varCount, withComponent, in);
  const bool componentsAscendingDesired = order.componentsAscendingDesired();
  const bool schreyering = order.schreyering();
  PolyRing ring(std::move(baseField), Monoid(std::move(order)));

  Processor processor(ring.monoid(), componentsAscendingDesired, schreyering);

  return std::make_pair(std::move(ring), std::move(processor));
}

void MathicIO::writeRing(
  const PolyRing& ring,
  const Processor& processor,
  const bool withComponent,
  std::ostream& out
){
  writeBaseField(ring.field(), out);
  out << ' ' << ring.varCount() << '\n';

  auto&& order = ring.monoid().makeOrder(
    processor.componentsAscendingDesired(),
    processor.schreyering()
  );
  writeOrder(order, withComponent, out);
}

auto MathicIO::readOrder(
  const VarIndex varCount,
  const bool withComponent,
  Scanner& in
) -> Order {
  const bool schreyering = in.match("schreyer");
  bool lexBaseOrder = !in.match("revlex") && in.match("lex");

  const auto gradingCount = in.readInteger<VarIndex>();
  bool componentsAscendingDesired = true;
  auto componentCompareIndex = Order::ComponentAfterBaseOrder;
  Gradings gradings(static_cast<size_t>(varCount) * gradingCount);
  size_t index = 0;
  for (VarIndex grading = 0; grading <  gradingCount; ++grading) {
    const bool com = withComponent && in.match("component");
    if (withComponent && (com || in.match("revcomponent"))) {
      if (!Monoid::HasComponent)
        in.reportError("Cannot specify component comparison for non-modules.");
      if (componentCompareIndex != Order::ComponentAfterBaseOrder)
        in.reportError("Component comparison must be specified at most once.");
      componentsAscendingDesired = com;
      componentCompareIndex = grading;
      index += varCount;
    } else {
      for (VarIndex i = 0; i < varCount; ++i, ++index)
        gradings[index] = in.readInteger<Exponent>();
    }
  }
  MATHICGB_ASSERT(index == gradings.size());

  const bool moreLex = in.match("_lex");
  if (moreLex || in.match("_revlex")) {
    lexBaseOrder = moreLex;
    const bool moreCom = withComponent && in.match("component");
    if (withComponent && (moreCom || in.match("revcomponent")))
      componentsAscendingDesired = moreCom;
  }

  Order order(
    varCount,
    std::move(gradings),
    lexBaseOrder ? Order::LexBaseOrder : Order::RevLexBaseOrder,
    componentCompareIndex,
    componentsAscendingDesired,
    schreyering
  );
  return std::move(order);
}

void MathicIO::writeOrder(
  const Order& order,
  const bool withComponent,
  std::ostream& out
) {
  MATHICGB_ASSERT(Monoid::HasComponent || !withComponent);

  const auto baseOrder =
    order.baseOrder() == Order::LexBaseOrder ? "lex" : "revlex";
  const auto componentOrder =
    order.componentsAscendingDesired() ? "component" : "revcomponent";

  if (order.schreyering())
    out << "schreyer ";
  out << baseOrder << ' ' << order.gradingCount() << '\n';
  for (VarIndex grading = 0; grading < order.gradingCount(); ++grading) {
    if (withComponent && grading == order.componentGradingIndex())
      out << ' ' << componentOrder;
    else {
      for (VarIndex var = 0; var < order.varCount(); ++var) {
        const auto index = var + grading * order.varCount();
        out << ' ' << unchar(order.gradings()[index]);
      }
    }
    out << '\n';
  }
  if (
    withComponent &&
    !order.componentsAscendingDesired() &&
    order.componentGradingIndex() == Order::ComponentAfterBaseOrder
  ) {
    out << " _" << baseOrder << "\n " << componentOrder << '\n';
  }
}

/*
auto MathicIO::readBasis() -> BasisData {
  typedef PolyRing::Monoid Monoid;
  typedef Monoid::VarIndex VarIndex;
  typedef Monoid::Gradings Gradings;

  // ** Read prime field
  const auto charac = mIn.readInteger<coefficient>();

  // ** Read 
  const auto varCount = mIn.readInteger<VarIndex>();
  const bool doSchreyer = mIn.match("schreyer");

  bool lexBaseOrder = false;
  if (mIn.match("lex"))
    lexBaseOrder = true;
  else if (mIn.match("revlex"))
    lexBaseOrder = false;

  const auto gradingCount = mIn.readInteger<VarIndex>();
  bool componentsAscendingDesired = true;
  auto componentCompareIndex = Order::ComponentAfterBaseOrder;
  Gradings gradings(static_cast<size_t>(varCount) * gradingCount);
  size_t index = 0;
  for (VarIndex grading = 0; grading <  gradingCount; ++grading) {
    const bool com = mIn.match("component");
    const bool revcom = !com && mIn.match("revcomponent");
    if (com || revcom) {
      if (!HasComponent)
        mIn.reportError("Cannot specify component comparison for non-modules.");
      componentsAscendingDesired = com;
      if (componentCompareIndex != Order::ComponentAfterBaseOrder)
        mIn.reportError("Component comparison must be specified at most once.");
      componentCompareIndex = grading;
      index += varCount;
    } else {
      for (VarIndex i = 0; i < varCount; ++i, ++index)
        gradings[index] = mIn.readInteger<Exponent>();
    }
  }
  MATHICGB_ASSERT(index == gradings.size());

  const bool moreRevlex = mIn.match("revlex");
  const bool moreLex = !moreRevlex && mIn.match("lex");
  if (moreRevlex || moreLex) {
    lexBaseOrder = moreLex;

    const bool moreCom = mIn.match("component");
    const bool moreRevcom = !moreCom && mIn.match("revcomponent");
    if (com || revcom)
      componentsAscendingDesired = moreCom;
  }

  auto ring = make_unique<PolyRing>(
    charac,
    Order(
      varCount,
      std::move(gradings),
      lexBaseOrder ? Order::LexBaseOrder : Order::RevLexBaseOrder,
      componentCompareIndex
    )
  );
  auto basis = make_unique<Basis>(*ring);
  auto processor = make_unique<MonoProcessor<Monoid>>(ring->monoid());
  processor->setComponentsAscendingDesired(componentsAscendingDesired);

  // ** Read polynomials
  const auto polyCount = mIn.readInteger<size_t>();
  for (size_t i = 0; i < polyCount; ++i) {
    auto p = readPolynomial();
    p->sortTermsDescending();
    basis->insert(std::move(p));
  }

  if (doSchreyer) {
    Monoid::MonoVector schreyer(ring->monoid());
    for (size_t gen = 0; gen < basis->size(); ++gen)
      schreyer.push_back(basis->getPoly(gen)->getLeadMonomial());
    processor->setModuleAdjustments(std::move(schreyer));
  }

  return std::make_tuple(
    std::move(ring),
    std::move(basis),
    std::move(processor)
  );
}

std::unique_ptr<Poly> MathicIO::readPolynomial(const PolyRing& ring) {
  Poly p(ring);
  if (mIn.match('0'))
    return;

  auto mono = ring.monoid().alloc();
  while (true) {
    if (!isdigit(next) && !isalpha(next) && next != '<')
      break;


    // read coefficient
    const bool negate = !mIn.match('+') && mIn.match('-');
    coefficient coef = negate ? -1 : 1;
    mIn.matchReadInteger(coef, negate);

    // read monomial
    read
    
    ring.monoid().setToIdentity(mono);
    

    const size_t firstLocation = monoms.size();
    monoms.resize(firstLocation + R->maxMonomialSize());
    monomial m = &monoms[firstLocation];
    if (isalpha(next) || next == '<')
      R->monomialParse(i, m);
    else
      R->monomialSetIdentity(m); // have to do this to set hash value
    next = i.peek();
    if (next == '>')
      i.get();
  }
}
*/

void MathicIO::readTerm(
  const PolyRing& ring,
  const bool readComponent,
  Coefficient& coef,
  MonoRef mono,
  Scanner& in
) {
  // ** Read coefficient, if any.
  const auto& field = ring.field();
  const bool negate = !in.match('+') && in.match('-');
  if (in.peekDigit())
    coef = field.toElement(in.readInteger<RawCoefficient>(negate));
  else if (negate)
    coef = field.minusOne();
  else
    coef = field.one();

  // ** Read monomial
  auto& monoid = ring.monoid();
  if (field.isOne(coef) && !in.peekAlpha()) {
    // Detect the monomial 1.
    monoid.setIdentity(mono);
    if (readComponent)
      this->readComponent(monoid, mono, in);
  } else
    readMonomial(monoid, readComponent, mono, in);
}

void MathicIO::readMonomial(
  const Monoid& monoid,
  const bool readComponent,
  MonoRef mono,
  Scanner& in
) {
  MATHICGB_ASSERT(!readComponent || Monoid::HasComponent);

  monoid.setIdentity(mono);
  if (in.peek() == '1') {
    const auto e = in.readInteger<Exponent>();
    if (e != 1) {
      std::ostringstream err;
      err << "Expected monomial, but got " << e << " (did you mean 1?).";
      in.reportError(err.str());
    }
  } else {
    bool sawSome = false;
    while (true) {
      const auto letterCount = 'z' - 'a' + 1;
      const auto letter = in.peek();
      
      VarIndex var;
      if ('a' <= letter && letter <= 'z')
        var = letter - 'a';
      else if ('A' <= letter && letter <= 'Z')
        var = (letter - 'A') + letterCount;
      else if (sawSome)
        break;
      else {
        std::ostringstream err;
        err << "Expected letter while reading monomial, but got '"
          << static_cast<char>(letter) << "'.";
        in.reportError(err.str());
        return;
      }
      in.get(); // skip past letter
      
      MATHICGB_ASSERT(var < 2 * letterCount);
      if (var >= monoid.varCount()) {
        std::ostringstream err;
        err << "Saw the variable " << static_cast<char>(letter)
          << ", but the monoid only has "
          << monoid.varCount() << " variables.";
        in.reportError(err.str());
        return;
      }
      if (monoid.exponent(mono, var) > static_cast<Exponent>(0)) {
        std::ostringstream err;
        err << "Variable " << static_cast<char>(letter) <<
          " must not be written twice in one monomial.";
        in.reportError(err.str());
      }
      
      if (in.peekDigit())
        monoid.setExponent(var, in.readInteger<Exponent>(), mono);
      else
        monoid.setExponent(var, static_cast<Exponent>(1), mono);
      sawSome = true;
    }
  }

  if (readComponent)
    this->readComponent(monoid, mono, in);
}

void MathicIO::readComponent(
  const Monoid& monoid,
  MonoRef mono,
  Scanner& in
) {
  MATHICGB_ASSERT(Monoid::HasComponent);
  in.expect('<');
  monoid.setComponent(in.readInteger<Exponent>(), mono);
  in.expect('>');
}

/// Print a monomial with no coefficient.
void MathicIO::writeMonomial(
  const Monoid& monoid,
  const bool writeComponent,
  ConstMonoRef mono,
  std::ostream& out
) {
  const auto letterCount = 'z' - 'a' + 1;

  bool printedSome = false;
  for (VarIndex var = 0; var < monoid.varCount(); ++var) {
    if (monoid.exponent(mono, var) == 0)
      continue;
    char letter;
    if (var < letterCount)
      letter = 'a' + static_cast<char>(var);
    else if (var < 2 * letterCount)
      letter = 'A' + (static_cast<char>(var) - letterCount);
    else {
      mathic::reportError("Too few letters in alphabet to print variable.");
      return;
    }
    printedSome = true;
    out << letter;
    if (monoid.exponent(mono, var) != 1)
      out << unchar(monoid.exponent(mono, var));
  }
  if (!printedSome)
    out << '1';
  if (writeComponent) {
    MATHICGB_ASSERT(Monoid::HasComponent);
    out << '<' << unchar(monoid.component(mono)) << '>';
  }
}

#endif
