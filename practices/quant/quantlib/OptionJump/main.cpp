/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2005, 2006, 2007, 2009 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

// #include <ql/qldefines.hpp>
// #if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
// #  include <ql/auto_link.hpp>
// #endif
// #include <ql/instruments/vanillaoption.hpp>
// #include <ql/pricingengines/vanilla/binomialengine.hpp>
// #include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
// #include <ql/pricingengines/vanilla/analytichestonengine.hpp>
// #include <ql/pricingengines/vanilla/baroneadesiwhaleyengine.hpp>
// #include <ql/pricingengines/vanilla/bjerksundstenslandengine.hpp>
// #include <ql/pricingengines/vanilla/batesengine.hpp>
// #include <ql/pricingengines/vanilla/integralengine.hpp>
// #include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
// #include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
// #include <ql/pricingengines/vanilla/mcamericanengine.hpp>
// #include <ql/pricingengines/vanilla/analyticeuropeanvasicekengine.hpp>
// #include <ql/time/calendars/target.hpp>
// #include <ql/utilities/dataformatters.hpp>
// #include <ql/models/shortrate/onefactormodels/vasicek.hpp>

#include <ql/quantlib.hpp>
#include <ql/processes/blackscholesprocess.hpp>

#include "business248.hpp"

#include <iostream>
#include <iomanip>

using namespace QuantLib;
// namespace ql = QuantLib;

class BSMJProcess : public BlackScholesMertonProcess {
public:
    BSMJProcess(
        const Handle<Quote>& x0,
        const Handle<YieldTermStructure>& dividendTS,
        const Handle<YieldTermStructure>& riskFreeTS,
        const Handle<BlackVolTermStructure>& blackVolTS,
        // Jump parameters
        const Real jumpint,
        const Real logJMean,
        const Real logJVol,
        const Real kappa, 
        // Arguments with default values
        const MersenneTwisterUniformRng& URng_Pt = MersenneTwisterUniformRng(SeedGenerator::instance().get()),
        const MersenneTwisterUniformRng& URng_Ji = MersenneTwisterUniformRng(SeedGenerator::instance().get()), // Jump size distribution
        const ext::shared_ptr<discretization>& d =
                ext::shared_ptr<discretization>(new EulerDiscretization)
        ):BlackScholesMertonProcess(x0, dividendTS, riskFreeTS, blackVolTS, d), jumpIntensity_(jumpint), logMeanJump_(logJMean), logJumpVolatility_(logJVol),
        kappa_(kappa), URng_Pt_(URng_Pt), URng_Ji_(URng_Ji){}

    // For simulation, we need to implement evolve(). Accordingly, the drift and diffusion terms should be adjusted.

    // Implement the adjustment for the drift term due to jump
    Real drift_adj(Time dt) const{
        return jumpIntensity_ * kappa_ * dt;
    }
    // Implement the adjustment for the diffusion term due to jump
    Real diffusion_adj(Time dt) const{
        InverseCumulativePoisson invP(jumpIntensity_ * dt); // It works
        InverseCumulativeNormal invN(logMeanJump_, logJumpVolatility_);
        Real total_jump = 0, Ji;

        Size Pt = static_cast<Size>(invP(URng_Pt_.next().value)); // Non-negative integer
        for (auto i=1; i<= Pt; i++){
            Ji = invN(URng_Ji_.next().value);
            total_jump += Ji;
        }
        return total_jump;
    }
    // Override the original evolve function by the modified evolve function with jump terms
    Real evolve(Time t0, Real x0, Time dt, Real dw) const {
        return apply(x0, discretization_->drift(*this, t0, x0, dt) - drift_adj(dt) +
            stdDeviation(t0, x0, dt)*dw + diffusion_adj(dt));
    }

private:
// Additional parameters for the jump diffusion modelling
    double jumpIntensity_;
    double logMeanJump_;
    double logJumpVolatility_;
    double kappa_; // Average jump size – 1
    // Uniform random number generators for the Monte Carlo Sim
    MersenneTwisterUniformRng URng_Pt_; // Poisson process for counting (lambda*dt)
    MersenneTwisterUniformRng URng_Ji_; // Jump size distribution

};

int main(int, char* []) {
    // set up dates
    Calendar calendar = TARGET();
    Date todaysDate(15, May, 1998);
    Date settlementDate(17, May, 1998);
    Settings::instance().evaluationDate() = todaysDate;

    // our options
    Option::Type type(Option::Put);
    Real underlying = 36;
    Real strike = 40;
    Spread dividendYield = 0.00;
    Rate riskFreeRate = 0.06;
    Volatility volatility = 0.20;
    Date maturity(17, May, 1999);
    DayCounter dayCounter = Actual365Fixed();

    std::cout << "Option type = "  << type << std::endl;
    std::cout << "Maturity = "        << maturity << std::endl;
    std::cout << "Underlying price = "        << underlying << std::endl;
    std::cout << "Strike = "                  << strike << std::endl;
    std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)
                << std::endl;
    std::cout << "Dividend yield = " << io::rate(dividendYield)
                << std::endl;
    std::cout << "Volatility = " << io::volatility(volatility)
                << std::endl;
    std::cout << std::endl;
    std::string method;
    std::cout << std::endl ;

    Handle<Quote> underlyingH(
        ext::shared_ptr<Quote>(new SimpleQuote(underlying)));

    // bootstrap the yield/dividend/vol curves
    Handle<YieldTermStructure> flatTermStructure(
        ext::shared_ptr<YieldTermStructure>(
            new FlatForward(settlementDate, riskFreeRate, dayCounter)));
    Handle<YieldTermStructure> flatDividendTS(
        ext::shared_ptr<YieldTermStructure>(
            new FlatForward(settlementDate, dividendYield, dayCounter)));
    Handle<BlackVolTermStructure> flatVolTS(
        ext::shared_ptr<BlackVolTermStructure>(
            new BlackConstantVol(settlementDate, calendar, volatility,
                                    dayCounter)));
    ext::shared_ptr<StrikedTypePayoff> payoff(
                                    new PlainVanillaPayoff(type, strike));

    Real jumpIntensity = 0.4736842;
    Real jumpVolatility = 0.4647414;
    Real jumpMean = log(0.9989906) - 0.5* jumpVolatility* jumpVolatility;
    Real kappa = 0.9989906 - 1; // The average jump size - 1

    // Construct the stock movement process object and
    // the pointer of the object

    boost::shared_ptr<BSMJProcess> bsmjProcess(new BSMJProcess(
        underlyingH, flatDividendTS, flatTermStructure, flatVolTS,
        jumpIntensity, jumpMean, jumpVolatility, kappa)
    );

    ext::shared_ptr<Exercise> europeanExercise(
                                    new EuropeanExercise(maturity));
    VanillaOption europeanOption(payoff, europeanExercise);

    Size timeSteps = 1;
    Size mcSeed = 42;    
    europeanOption.setPricingEngine(MakeMCEuropeanEngine<PseudoRandom>(bsmjProcess)
        .withSteps(timeSteps)
        .withAbsoluteTolerance(0.02)
        .withSeed(mcSeed));
    return 0;

}
