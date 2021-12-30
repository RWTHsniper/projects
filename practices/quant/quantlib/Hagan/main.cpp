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

#include "main.hpp"

namespace ql = QuantLib;
namespace pt = boost::property_tree;
using namespace QuantLib;

int main(int, char* []) {

    // Read input file
    std::string file_dir(__FILE__);
    std::string source_dir(file_dir);
    eraseSubStr(source_dir, "main.cpp"); // remove main.cpp from the file's path
    std::string data_dir(source_dir+"data/");

    // Short alias for this namespace
    // Create a root
    pt::ptree KRWIRS;
    // Load the json file in this ptree
    pt::read_json(data_dir+"KRWIRS.json", KRWIRS);
    const size_t& swap_size = KRWIRS.get("size", 0); // Length of the term-structure

    // Set up dates
    Date todaysDate(Date(DateParser::parseFormatted(KRWIRS.get("base_date", "NONE"),"%Y-%m-%d")));
    Settings::instance().evaluationDate() = todaysDate;
    std::cout << "\nThe evaluation date for this example is set to "
                << Settings::instance().evaluationDate() << std::endl;
    std::vector<Date> tenorDates;
    tenorDates.reserve(swap_size);
    for (const pt::ptree::value_type &tenor_date : KRWIRS.get_child("tenor_dates")){
        tenorDates.emplace_back(Date(DateParser::parseFormatted(tenor_date.second.data(),"%Y-%m-%d")));
    }
    std::vector<Period> tenors;
    tenors.reserve(swap_size);
    for (const pt::ptree::value_type &tenor : KRWIRS.get_child("tenors")){
        tenors.emplace_back(Period(PeriodParser::parse(tenor.second.data())));
    }
    std::vector<double> rates;
    rates.reserve(swap_size);
    for (const pt::ptree::value_type &rate : KRWIRS.get_child("rates")){
        rates.emplace_back(stod(rate.second.data()));
    }

    // Term-structure in QuantLib
    std::vector<ext::shared_ptr<RateHelper>> swapHelpers;
    swapHelpers.reserve(swap_size);
    for (size_t i=0; i<tenors.size(); i++){
        // I tried EuriborSwapIsdaFixA, GbpLiborSwapIsdaFix, and JpyLiborSwapIsdaFixAm
        // The maximum difference is 13.26 bp. Thus, I just choose EuriborSwapIsdaFixA tentatively.
        ext::shared_ptr<SwapIndex> swapInd(new EuriborSwapIsdaFixA(tenors[i])); // tentatively chosen.
        ext::shared_ptr<RateHelper> helper(new SwapRateHelper(rates[i],swapInd));
        swapHelpers.emplace_back(helper);
    }
    // Discount: the rate which should be constructed
    // Possible interpolators: Cubic, LogLinear, Linear
    double tolerance = 1.0e-15;
    // ext::shared_ptr<YieldTermStructure> yieldCurve(new PiecewiseYieldCurve<Discount, Cubic>(todaysDate, swapHelpers, Actual365Fixed(),PiecewiseYieldCurve<Discount, Cubic>::bootstrap_type(tolerance)));
    // ext::shared_ptr<YieldTermStructure> yieldCurve(new PiecewiseYieldCurve<Discount, LogLinear>(todaysDate, swapHelpers, Actual365Fixed(),PiecewiseYieldCurve<Discount, LogLinear>::bootstrap_type(tolerance)));
    ext::shared_ptr<YieldTermStructure> yieldCurve(new PiecewiseYieldCurve<Discount, Linear>(todaysDate, swapHelpers, Actual365Fixed(),PiecewiseYieldCurve<Discount, Linear>::bootstrap_type(tolerance)));
    std::cout << yieldCurve->discount(tenorDates[5]) << std::endl;
    for (size_t i=0; i<tenorDates.size(); i++){
        const auto& date = tenorDates[i];
        std::cout << date << " " << yieldCurve->discount(date) << std::endl;
    }

    // Black's normal volatility model in QuantLib
    

    // SwapRateHelper()
    /* 
    Cpp

             PiecewiseYieldCurve(
             const Date& referenceDate,
             std::vector<ext::shared_ptr<typename Traits::helper> > instruments,
             const DayCounter& dayCounter,
             const std::vector<Handle<Quote> >& jumps = std::vector<Handle<Quote> >(),
             const std::vector<Date>& jumpDates = std::vector<Date>(),
             const Interpolator& i = Interpolator(),
             bootstrap_type bootstrap = bootstrap_type())
         : base_curve(referenceDate, dayCounter, jumps, jumpDates, i),
           instruments_(std::move(instruments)), accuracy_(1.0e-12),
           bootstrap_(std::move(bootstrap)) {
             bootstrap_.setup(this);
         }

    Python code
    ql.SwapRateHelper(0.06, ql.EuriborSwapIsdaFixA(ql.Period('1y')))
    curve = ql.PiecewiseLogLinearDiscount(ql.Date(15,6,2020), helpers, ql.Actual360())



SwapRateHelper	(	Rate 	rate,
const ext::shared_ptr< SwapIndex > & 	swapIndex,
Handle< Quote > 	spread = Handle<Quote>(),
const Period & 	fwdStart = 0 * Days,
Handle< YieldTermStructure > 	discountingCurve = Handle<YieldTermStructure>(),
Pillar::Choice 	pillar = Pillar::LastRelevantDate,
Date 	customPillarDate = Date(),
bool 	endOfMonth = false 
)	
    */
    // PiecewiseYieldCurve (const Date &referenceDate, std::vector< ext::shared_ptr< typename Traits::helper > > instruments, const DayCounter &dayCounter, bootstrap_type bootstrap)



    // // Create a root
    // pt::ptree swaption_vol;
    // // Load the json file in this ptree
    // pt::read_json(data_dir+"swaption_vol.json", swaption_vol);
    // expiry: option. tenor: swap


    /*
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
*/
    return 0;

}
