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
    const Size& swap_size = KRWIRS.get("size", 0); // Length of the term-structure

    // Set up dates
    DayCounter actualDC = Actual365Fixed();
    Date todaysDate(Date(DateParser::parseFormatted(KRWIRS.get("base_date", "NONE"),"%Y-%m-%d")));
    Settings::instance().evaluationDate() = todaysDate;
    std::cout << "\nThe evaluation date for this example is set to "
                << Settings::instance().evaluationDate() << std::endl;
    std::vector<Date> tenorDates;
    tenorDates.reserve(swap_size);
    for (const pt::ptree::value_type &tenor_date : KRWIRS.get_child("tenor_dates")){
        tenorDates.emplace_back(DateParser::parseFormatted(tenor_date.second.data(),"%Y-%m-%d"));
    }
    std::vector<Period> tenors;
    tenors.reserve(swap_size);
    for (const pt::ptree::value_type &tenor : KRWIRS.get_child("tenors")){
        tenors.emplace_back(PeriodParser::parse(tenor.second.data()));
    }
    std::vector<double> rates;
    rates.reserve(swap_size);
    for (const pt::ptree::value_type &rate : KRWIRS.get_child("rates")){
        rates.emplace_back(stod(rate.second.data()));
    }

    // Term-structure in QuantLib
    Calendar korCal = SouthKorea(); 
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
    ext::shared_ptr<YieldTermStructure> yieldCurve(new PiecewiseYieldCurve<Discount, LogLinear>(todaysDate, swapHelpers, Actual365Fixed(),PiecewiseYieldCurve<Discount, LogLinear>::bootstrap_type(tolerance)));
    // ext::shared_ptr<YieldTermStructure> yieldCurve(new PiecewiseYieldCurve<Discount, Linear>(todaysDate, swapHelpers, Actual365Fixed(),PiecewiseYieldCurve<Discount, Linear>::bootstrap_type(tolerance)));
    for (size_t i=0; i<tenorDates.size(); i++){
        const auto& date = tenorDates[i];
        std::cout << date << " " << yieldCurve->discount(date) << std::endl; // DiscountFactor
        std::cout << "fwd " << Time(date-todaysDate)  << " " << yieldCurve->forwardRate(todaysDate, date, Actual365Fixed(), Continuous) << std::endl; // InterestRate
    }

    // Bachilier's normal volatility model in QuantLib
    Swap::Type swapType = Swap::Payer; // type of swaption
    ext::shared_ptr<Quote> vol(new SimpleQuote(0.2)); // 20% of volatile for a test
    InterestRate fwdRate = yieldCurve->forwardRate(0.0,6.0,Simple); // compute forward rate from the yield curve
    ext::shared_ptr<Quote> flatRate(new SimpleQuote(fwdRate.rate())); // discounting
    Handle<YieldTermStructure> rhTermStructure(
    ext::make_shared<FlatForward>(
                todaysDate, Handle<Quote>(flatRate), // assume todaysDate = settlementDate
                                Actual365Fixed()));
                                ext::shared_ptr<IborIndex> indexThreeMonths(new Euribor3M(rhTermStructure));
    ext::shared_ptr<BlackCalibrationHelper> swaption(new SwaptionHelper(
                                                        Period(6, Years), // maturity Period
                                                        Period(1, Years), // length Period
                                                        Handle<Quote>(vol), // 
                                                        indexThreeMonths,
                                                        indexThreeMonths->tenor(),
                                                        Actual365Fixed(), // Actual365Fixed() is more exact than indexThreeMonths->dayCounter()
                                                        Actual365Fixed(), // Actual365Fixed()
                                                        rhTermStructure,
                                                        BlackCalibrationHelper::RelativePriceError,
                                                        Null<Real>(),
                                                        1.0,
                                                        Normal));
    ext::shared_ptr<PricingEngine> bachelierEngine(new BachelierSwaptionEngine(rhTermStructure, Handle<Quote>(ext::shared_ptr<Quote>(vol))));
    swaption->setPricingEngine(bachelierEngine);

    Real npv = swaption->modelValue();
    std::cout << "npv: " << npv << std::endl;
    Volatility implied = swaption->impliedVolatility(npv, 1e-4,1000, 0.05, 0.50);
    std::cout << "ivol: " << implied << std::endl;

    // Create a root
    pt::ptree swaption_vol;
    // Load the json file in this ptree
    pt::read_json(data_dir+"swaption_vol.json", swaption_vol); // expiry: option. tenor: swap
    const Size& expiry_size = swaption_vol.get("expiry_size", 0); // Length of a swaption's expiry
    const Size& tenor_size = swaption_vol.get("tenor_size", 0); // Length of a swaption's tenor
    std::cout << "expiry and tenor sizes " << expiry_size << " " << tenor_size << std::endl;
    std::shared_ptr<std::vector<Period>> swaptionExpiry = std::make_shared<std::vector<Period>>();
    swaptionExpiry->reserve(expiry_size);
    for (const pt::ptree::value_type &expiry : swaption_vol.get_child("expiry")){
        swaptionExpiry->emplace_back(stoi(expiry.second.data()), Years);
    }
    std::shared_ptr<std::vector<Period>> swaptionTenor = std::make_shared<std::vector<Period>>();
    swaptionTenor->reserve(tenor_size);
    for (const pt::ptree::value_type &tenor : swaption_vol.get_child("tenor")){
        swaptionTenor->emplace_back(stoi(tenor.second.data()), Years);
    }
    // read swaption ATM volatility matrix
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> swaptionVolMat(expiry_size, tenor_size); // row-major matrix with dynamic allocation
    std::shared_ptr<Eigen::MatrixXd> swaptionVolMat = std::make_shared<Eigen::MatrixXd>(expiry_size, tenor_size);
    // Eigen::MatrixXd swaptionVolMat(expiry_size, tenor_size); // row-major matrix with dynamic allocation
    std::array<size_t, 2> counter{0,0};
    for (const pt::ptree::value_type &quote : swaption_vol.get_child("quote")){
        (*swaptionVolMat)(counter[0], counter[1]) = stod(quote.second.data());
        counter[1]++; 
        if (counter[1] % tenor_size == 0){
            counter[0]++;
            counter[1] = 0;
        }
    }

    std::cout << *swaptionVolMat << std::endl;

    // Read hagan's LGM model
    pt::ptree lgm;
    pt::read_json(data_dir+"LGM.json", lgm); // expiry: option. tenor: swap
    const size_t& nFactor = lgm.get("numFactors", 0); // number of factors
    counter[0] = 0; counter[1] = 0; // initialize the counter
    Eigen::MatrixXd corrMat(nFactor, nFactor);
    for (const pt::ptree::value_type &corr : lgm.get_child("corrUpper")){ // Upper triangular part in the correlation matrix
        corrMat(counter[0], counter[1]) = stod(corr.second.data());
        corrMat(counter[1], counter[0]) = corrMat(counter[0], counter[1]); // symmetric
        counter[1]++; 
        if (counter[1] % nFactor == 0){
            counter[0]++; // next row
            counter[1] = counter[0]; // diagonal
        }
    }

    std::cout << "This is corrMat " << std::endl << corrMat << std::endl;
    // corrMat << 1.0,0.5,0.5,1.0;
    StochasticModel::HaganNF haganModel(yieldCurve, nFactor, corrMat);

    // testModel(); // test classes in model.hpp

    // Eigen::VectorXd a(3); a<< 1,2,3;
    // auto b = Eigen::square(a); std::cout << "a " << b <<std::endl;

    // Simulation
    // const size_t numPaths = 10;
    const double T = 2;
    const size_t numPaths = 100; // test
    // const size_t numSteps = 365*T; // daily
    const size_t numSteps = 52*T; // weekly

    // Initialization for Brownian motion
    std::vector<Eigen::MatrixXd> dWIndep; // [numSteps](nFactor, numPaths)
    buildBrownianMotion(dWIndep, numSteps, nFactor, numPaths);
    const double dt = T/(numSteps);
    const double t_i = 0.0;
    double t = t;
    // std::cout << dWIndep[0] << std::endl;
    // std::cout << dWIndep[numSteps-1] << std::endl;
    // Eigen::VectorXd x(nFactor); x.setZero();
    std::vector<Eigen::MatrixXd> x;
    buildVectofMat(x, numSteps+1, nFactor, numPaths); // 0th element in vector represents the initial values
    Eigen::MatrixXd dw(nFactor, numPaths);
    // Eigen::VectorXd dw(nFactor);
    const Eigen::MatrixXd& lowerMat_ = haganModel.getLowerMat();
    // std::cout << "lower mat " << lowerMat_ << std::endl;
    for (size_t i=0; i < numSteps ;i++){
        // dw = lowerMat_ * dWIndep[i];
        dw = lowerMat_; dw *= dWIndep[i];
        haganModel.evolve(x[i+1], t, x[i], dt, dw);
        t+= dt;
        if (i==numSteps-1) saveData(source_dir+"output/dw.csv", dw);
    }

    saveData(source_dir+"output/matrix.csv", x[numSteps]); // save matrix in output folder
    std::shared_ptr<Eigen::MatrixXd> r =  haganModel.computeInterestRate(t_i, dt, numPaths, numSteps, x); // size of (numPaths, numSteps+1)
    saveData(source_dir+"output/r.csv", (*r)); // save interest rate paths

    ql::Period today(0, ql::Years); // Current time for evaluation
    double iVol = haganModel.impliedVol(today, (*swaptionExpiry)[0], (*swaptionTenor)[0], 0.25, ql::Normal);
    std::cout << "Ivol from HaganModel " << iVol << std::endl;


    Optimizer::testBoothFun();

    /* Next steps
        - Think about pricing swaptions. (analytic)
        - Think about pricing swaptions. (MC)
    */


        // x[0] = spotQuote -> value ();
        // x[1] = v0;
        // ql::Size numVals =10;
        // for ( ql::Size j=1;j <= numVals ;++j){
        //     dw [0]= bmGauss . next (). value ;
        //     dw [1]= bmGauss . next (). value ;
        //     x= hestonProcess -> evolve (t,x,dt ,dw );
        //     std :: cout << " Time : " << t+dt << ", S_t : " << x[0] << ", V_t : " << x [1] << std :: endl ;
        //     t+= dt;
        // }


/*
I worked on defining it, but not going to need it
    std::vector<std::vector<Handle<Quote>>> swaptionQuote;
    swaptionQuote.reserve(expiry_size);
    std::vector<Handle<Quote>> rowVec;
    rowVec.reserve(tenor_size);
    for (const pt::ptree::value_type &quote : swaption_vol.get_child("quote")){
        // std::cout << stod(quote.second.data()) << std::endl;
        rowVec.emplace_back(ext::shared_ptr<Quote>(new SimpleQuote(stod(quote.second.data()))));
        if (rowVec.size() == tenor_size){
            swaptionQuote.emplace_back(rowVec);
            rowVec.clear();
            rowVec.reserve(tenor_size);
        }
    }
    auto iVolMat = SwaptionVolatilityMatrix(todaysDate,
                                            korCal, 
                                            ModifiedFollowing, 
                                            swaptionExpiry,
                                            swaptionTenor,
                                            swaptionQuote,
                                            actualDC,
                                            false,
                                            Normal);
    std::cout << iVolMat[0][0] << std::endl;
*/
    // std::cout << iVolMat << " iVolMat" << std::endl;

/*
mySwaption test
    // Swap::Type swapType = Swap::Receiver;
    // Make a discount factor vector
    std::vector<Real> discountFactors;
    discountFactors.reserve(4);
    discountFactors.emplace_back(yieldCurve->discount(6.25));
    discountFactors.emplace_back(yieldCurve->discount(6.5));
    discountFactors.emplace_back(yieldCurve->discount(6.75));
    discountFactors.emplace_back(yieldCurve->discount(7.0));
    std::cout << "discount factors " << std::endl;
    for (const auto& df : discountFactors){
        std::cout << "df " << df << std::endl;
    }
    std::cout << "df at 6: " << yieldCurve->discount(6.0) << std::endl;
    std::vector<Real> prices;
    prices.reserve(5);
    prices.emplace_back(bachelierATMSwaption(0.0, 6, 0.25, 0.0, swapType, discountFactors));
    prices.emplace_back(bachelierATMSwaption(0.0, 6, 0.25, 0.05, swapType, discountFactors));
    prices.emplace_back(bachelierATMSwaption(0.0, 6, 0.25, 0.10, swapType, discountFactors));
    prices.emplace_back(bachelierATMSwaption(0.0, 6, 0.25, 0.15, swapType, discountFactors));
    prices.emplace_back(bachelierATMSwaption(0.0, 6, 0.25, 0.20, swapType, discountFactors));
    for (const auto& p : prices){
        std::cout << "p " << p << std::endl;
    }
    */

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
