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
// using namespace QuantLib;

int main(int, char* []) {

    // Read input file
    std::string file_dir(__FILE__);
    std::string source_dir(file_dir);
    eraseSubStr(source_dir, "main.cpp"); // remove main.cpp from the file's path
    std::string data_dir(source_dir+"data/");
    std::string output_dir(source_dir+"output/");

    // Short alias for this namespace
    // Create a root
    pt::ptree KRWIRS;
    // Load the json file in this ptree
    pt::read_json(data_dir+"KRWIRS.json", KRWIRS);
    const ql::Size& swap_size = KRWIRS.get("size", 0); // Length of the term-structure

    // Set up dates
    ql::DayCounter actualDC = ql::Actual365Fixed();
    ql::Date todaysDate(ql::Date(ql::DateParser::parseFormatted(KRWIRS.get("base_date", "NONE"),"%Y-%m-%d")));
    ql::Settings::instance().evaluationDate() = todaysDate;
    std::cout << "\nThe evaluation date for this example is set to "
                << ql::Settings::instance().evaluationDate() << std::endl;
    std::vector<ql::Date> tenorDates;
    tenorDates.reserve(swap_size);
    for (const pt::ptree::value_type &tenor_date : KRWIRS.get_child("tenor_dates")){
        tenorDates.emplace_back(ql::DateParser::parseFormatted(tenor_date.second.data(),"%Y-%m-%d"));
    }
    std::vector<ql::Period> tenors;
    tenors.reserve(swap_size);
    for (const pt::ptree::value_type &tenor : KRWIRS.get_child("tenors")){
        tenors.emplace_back(ql::PeriodParser::parse(tenor.second.data()));
    }
    std::vector<double> rates;
    rates.reserve(swap_size);
    for (const pt::ptree::value_type &rate : KRWIRS.get_child("rates")){
        rates.emplace_back(stod(rate.second.data()));
    }

    // Term-structure in QuantLib
    ql::Calendar korCal = ql::SouthKorea(); 
    std::vector<boost::shared_ptr<ql::RateHelper>> swapHelpers;
    swapHelpers.reserve(swap_size);
    for (size_t i=0; i<tenors.size(); i++){
        // I tried EuriborSwapIsdaFixA, GbpLiborSwapIsdaFix, and JpyLiborSwapIsdaFixAm
        // The maximum difference is 13.26 bp. Thus, I just choose EuriborSwapIsdaFixA tentatively.
        boost::shared_ptr<SwapIndex> swapInd(new EuriborSwapIsdaFixA(tenors[i])); // tentatively chosen.
        boost::shared_ptr<RateHelper> helper(new SwapRateHelper(rates[i],swapInd));
        swapHelpers.emplace_back(helper);
    }
    // Discount: the rate which should be constructed
    // Possible interpolators: Cubic, LogLinear, Linear
    double tolerance = 1.0e-15;
    // ext::shared_ptr<YieldTermStructure> yieldCurve(new PiecewiseYieldCurve<Discount, Cubic>(todaysDate, swapHelpers, Actual365Fixed(),PiecewiseYieldCurve<Discount, Cubic>::bootstrap_type(tolerance)));
    boost::shared_ptr<ql::YieldTermStructure> yieldCurve(new ql::PiecewiseYieldCurve<ql::Discount, ql::LogLinear>(
                            todaysDate, swapHelpers, ql::Actual365Fixed(),ql::PiecewiseYieldCurve<ql::Discount, ql::LogLinear>::bootstrap_type(tolerance)));
    // ext::shared_ptr<YieldTermStructure> yieldCurve(new PiecewiseYieldCurve<Discount, Linear>(todaysDate, swapHelpers, Actual365Fixed(),PiecewiseYieldCurve<Discount, Linear>::bootstrap_type(tolerance)));
    // for (size_t i=0; i<tenorDates.size(); i++){
    //     const auto& date = tenorDates[i];
    //     std::cout << date << " " << yieldCurve->discount(date) << std::endl; // DiscountFactor
    //     std::cout << "fwd " << Time(date-todaysDate)  << " " << yieldCurve->forwardRate(todaysDate, date, Actual365Fixed(), Continuous) << std::endl; // InterestRate
    // }

    // Bachilier's normal volatility model in QuantLib
    // ql::Swap::Type swapType = ql::Swap::Payer; // type of swaption
    // boost::shared_ptr<ql::Quote> vol(new ql::SimpleQuote(0.2)); // 20% of volatile for a test
    // ql::InterestRate fwdRate = yieldCurve->forwardRate(0.0,6.0,Simple); // compute forward rate from the yield curve
    // boost::shared_ptr<ql::Quote> flatRate(new ql::SimpleQuote(fwdRate.rate())); // discounting
    // ql::Handle<ql::YieldTermStructure> rhTermStructure(
    // boost::make_shared<ql::FlatForward>(
    //             todaysDate, ql::Handle<ql::Quote>(flatRate), // assume todaysDate = settlementDate
    //                             ql::Actual365Fixed()));
    //                             boost::shared_ptr<IborIndex> indexThreeMonths(new ql::Euribor3M(rhTermStructure));
    // boost::shared_ptr<BlackCalibrationHelper> swaption(new ql::SwaptionHelper(
    //                                                     ql::Period(6, ql::Years), // maturity Period
    //                                                     ql::Period(1, ql::Years), // length Period
    //                                                     ql::Handle<ql::Quote>(vol), // 
    //                                                     indexThreeMonths,
    //                                                     indexThreeMonths->tenor(),
    //                                                     ql::Actual365Fixed(), // Actual365Fixed() is more exact than indexThreeMonths->dayCounter()
    //                                                     ql::Actual365Fixed(), // Actual365Fixed()
    //                                                     rhTermStructure,
    //                                                     ql::BlackCalibrationHelper::RelativePriceError,
    //                                                     ql::Null<Real>(),
    //                                                     1.0,
    //                                                     ql::Normal));
    // boost::shared_ptr<ql::PricingEngine> bachelierEngine(new ql::BachelierSwaptionEngine(rhTermStructure, ql::Handle<ql::Quote>(boost::shared_ptr<ql::Quote>(vol))));
    // swaption->setPricingEngine(bachelierEngine);

    // ql::Real npv = swaption->modelValue();
    // std::cout << "npv: " << npv << std::endl;
    // ql::Volatility implied = swaption->impliedVolatility(npv, 1e-4,1000, 0.05, 0.50);
    // std::cout << "ivol: " << implied << std::endl;

    // Eigen::VectorXd discountFactors(4);
    // discountFactors[0] = yieldCurve->discount(6.25);
    // discountFactors[1] = yieldCurve->discount(6.5);
    // discountFactors[2] = yieldCurve->discount(6.75);
    // discountFactors[3] = yieldCurve->discount(7.0);
    // double myIvol = StochasticModel::computeIvol(npv, 0.0, 6.0, 0.25, discountFactors, 0.0001);
    // std::cout << "My Ivol " << myIvol << std::endl;

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
    // row major like C
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
    const size_t& alpOrder = lgm.get("alpOrder", 0); // number of factors
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

    std::cout << "This is correlation matrix between factors " << std::endl << corrMat << std::endl;
    StochasticModel::HaganNF haganModel(yieldCurve, nFactor, corrMat, alpOrder);


    if(lgm.get("calibrate", false)){
        haganModel.calibrate(swaptionExpiry, swaptionTenor, swaptionVolMat); // skip calibration because I am working on validation
    }

    Eigen::MatrixXd computedIvol(swaptionExpiry->size(), swaptionTenor->size());
    for (size_t i=0; i<swaptionExpiry->size(); i++){
        for (size_t j=0; j<swaptionTenor->size(); j++){
            computedIvol(i,j) = haganModel.impliedVol((*swaptionExpiry)[i], (*swaptionTenor)[j], 0.25, ql::Normal);
        }
    }

    // Code for testing
    // double swapPeriod = 0.25;
    // Eigen::MatrixXd myIvols = haganModel.computeSwaptionMC(swaptionExpiry, swaptionTenor, x, swapPeriod, dt);
    // std::cout << "myIvol from sim" << std::endl << myIvols << std::endl;
    // std::cout << "computed" << std::endl << computedIvol << std::endl;
    // Write output as JSON
    std::string filePath = output_dir+"computedVol.json";
    writeOutputJSON(filePath, swaptionExpiry, swaptionTenor, computedIvol);
    std::cout << "computed Volatility " << std::endl << computedIvol << std::endl;
    std::cout << "Difference " << std::endl << computedIvol - *swaptionVolMat << std::endl;


    // Run simulation for pricing
    // const size_t numPaths = 10;
    const double T = lgm.get("T", 0);
    const size_t numPaths = lgm.get("numPaths", 0); // test
    const size_t numAnnualSteps = lgm.get("numAnnualSteps", 0); // 365: daily, 52: weekly, 12: monthly, 4: quarterly
    const size_t numSteps = numAnnualSteps*T; // monthly

    // Initialization for Brownian motion
    std::vector<Eigen::MatrixXd> dWIndep; // [numSteps](nFactor, numPaths)
    Simulation::buildBrownianMotion(dWIndep, numSteps, nFactor, numPaths);
    const double dt = T/(numSteps); // time-step size in the simulation
    const double t_i = 0.0; // start time in simulation
    double t = t_i;
    std::vector<Eigen::MatrixXd> x;
    Simulation::buildVectofMat(x, numSteps+1, nFactor, numPaths); // 0th element in vector represents the initial values
    Eigen::MatrixXd dw(nFactor, numPaths);
    // Eigen::VectorXd dw(nFactor);
    const Eigen::MatrixXd& lowerMat_ = haganModel.getLowerMat();
    // std::cout << "lower mat " << lowerMat_ << std::endl;
    for (size_t i=0; i < numSteps ;i++){
        // dw = lowerMat_ * dWIndep[i];
        dw = lowerMat_ * dWIndep[i];
        haganModel.evolve(x[i+1], t, x[i], dt, dw);
        t+= dt;
        // if (i==numSteps-1) Simulation::saveData(output_dir+"dw.csv", dw); // for a test, save Brownian motion
    }
    std::cout << "Running Monte-Carlo simulation is complete" << std::endl;

    // Load the json file in this ptree
    pt::ptree pricing;
    pt::read_json(data_dir+"pricing.json", pricing);
    const double& noteExpiry = pricing.get("expiry", 0.0); // number of factors
    const double& notecouponPeriod = pricing.get("couponPeriod", 0.0); // number of factors
    const double& notecouponRate = pricing.get("couponRate", 0.0); // number of factors
    const double& noteTau1 = pricing.get("tau1", 0.0); // Maturity for the first
    const double& noteFreq1 = pricing.get("freq1", 0.0); // coupon payment frequency for the first
    const double& noteTau2 = pricing.get("tau2", 0.0); // Maturity for the second
    const double& noteFreq2 = pricing.get("freq2", 0.0); // coupon payment frequency for the second
    // Compute spread_5_10
    const size_t& noteSteps = static_cast<size_t>(noteExpiry/dt);
    Eigen::MatrixXd spread_5_10(numPaths, noteSteps+1); spread_5_10.setZero();
    std::cout << "Start computing the spread between " << noteTau1 << "Y and " << noteTau2 << "Y IRS rates" << std::endl;
    haganModel.computeSpread(spread_5_10, x, dt, noteTau1, noteFreq1, noteTau2, noteFreq2, noteExpiry);
    // Compute the spread accrual...
    std::cout << "Start pricing a Spread Range Accrual Note " << std::endl;
    double SRAN_price = haganModel.computeSRAN(spread_5_10, x, notecouponRate, notecouponPeriod, dt);
    std::cout << "Spread Range Accrual Note's price" << std::endl << SRAN_price << std::endl;

    // Compute interest rate and save results
    std::shared_ptr<Eigen::MatrixXd> r =  haganModel.computeInterestRate(t_i, dt, numPaths, numSteps, x); // size of (numPaths, numSteps+1)
    std::shared_ptr<Eigen::MatrixXd> M =  Simulation::computeMSA(r, dt); // size of (numPaths, numSteps+1)
    Simulation::saveData(output_dir+"matrix.csv", x[numSteps]); // save matrix in output folder
    Simulation::saveData(output_dir+"r.csv", (*r)); // save interest rate paths
    Simulation::saveData(output_dir+"M.csv", (*M)); // save money-savings-account paths
    Simulation::saveData(output_dir+"spread_5_10.csv", spread_5_10); // save money-savings-account paths




    return 0;

}
