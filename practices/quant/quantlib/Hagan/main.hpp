#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <iostream>
#include <string> // stoi, stod
#include <typeinfo>

#include <Eigen/Dense> // Dense matrices
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <ql/quantlib.hpp>
#include <ql/time/calendar.hpp>
#include <ql/pricingengine.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>
#include <ql/utilities/dataparsers.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolmatrix.hpp>

#include "stochasticmodel.hpp"
#include "model.hpp"
#include "utils.hpp"

using namespace QuantLib;

/*
* Bachelier's Normal model for pricing swaptions
*/
Real bachelierATMSwaption(Real t, Real T0, Real tau, Real sig, Swap::Type type, std::vector<Real>& discountFactors, Real notional=1.0){
/*
Arguments
- T0: expiration/reset date 
- t: current time
- sig: implied volatility
- type: type of a swaption (Payer or Receiver)
- discountFactors: vector of discount factors at coupon payments
- notional: notional amount
*/
    Real payer_val{0.0};
    // Since the goal is to compute ATM swaption, R_swap(t) = K
    Real sqrt_t = sqrt(T0-t); // sqrt of time to expiry
    Real pdf_0{0.3989422804014327}; // pdf of Gaussian distribution at 0
    payer_val = std::accumulate(discountFactors.begin(), discountFactors.end(), 0.0);
    payer_val *= notional * tau * sig * sqrt_t * pdf_0;

    switch (type) {
        case Swap::Payer: return payer_val;
        case Swap::Receiver: return -payer_val;
        default: std::cout << "d"; exit(-1); // there are no applicable constant_expressions 
                                   // therefore default is executed
    }

}

void testModel(){
    double tol = 1e-8;
    Eigen::VectorXd xVals(5);
    Eigen::VectorXd yVals(5);
    xVals << 1,2,3,4,5;
    yVals << 1,4,9,16,25;
    Model::PolyFunc polyCurve(xVals, yVals, 2);
    for (size_t i=0; i<xVals.size(); i++){
        std::cout << polyCurve.evaluate(xVals[i]) << " " << yVals[i] << std::endl;
        assert(tol >= std::abs(polyCurve.evaluate(xVals[i]) - yVals[i])); // should be less than a tolerance
    }
    std::cout << "+1" << std::endl;
    polyCurve += 1;
    for (size_t i=0; i<xVals.size(); i++){
        std::cout << polyCurve.evaluate(xVals[i]) << " " << yVals[i] << std::endl;
        assert(tol >= std::abs(polyCurve.evaluate(xVals[i]) - (yVals[i]+1))); // should be less than a tolerance
    }
    std::cout << "Integration of x^2 + 1" << std::endl; // 1/3*x^3 + x 
    for (size_t i=0; i<xVals.size(); i++){
        std::cout << xVals[i] << " " << polyCurve.evalInt(0.0, xVals[i]) << std::endl;
        assert(tol >= std::abs(polyCurve.evalInt(0.0, xVals[i]) - (1/3.0 * pow(xVals[i],3) + xVals[i])));
    }
    Model::PolyFunc polyCurve2(polyCurve);
    std::cout << "add1 " << &(polyCurve) << std::endl;
    std::cout << "add2 " << &(polyCurve2) << std::endl;
    Model::PolyFunc polyCurve3 = polyCurve + polyCurve2; // (x^2+1)*2
    std::cout << "add3 " << &(polyCurve3) << std::endl;
    for (size_t i=0; i<xVals.size(); i++){
        std::cout << polyCurve3.evaluate(xVals[i]) << " " << yVals[i] << std::endl;
        assert(tol >= std::abs(polyCurve3.evaluate(xVals[i]) - 2.0*(pow(xVals[i],2)+1)));
    }
    Model::PolyFunc polyCurve4 = polyCurve*polyCurve2; // (x^2+1)^2
    for (size_t i=0; i<xVals.size(); i++){
        std::cout << polyCurve4.evaluate(xVals[i]) << " " << yVals[i] << std::endl;
        assert(tol >= std::abs(polyCurve4.evaluate(xVals[i]) - pow(pow(xVals[i],2)+1,2)));
    }

    Eigen::VectorXd params(3);
    params << 1.0, 1.0, 1.0;
    Eigen::VectorXd yExpVals(5);
    yExpVals << exp(1),exp(2),exp(3),exp(4),exp(5);
    Model::ExpFunc expCurve(xVals, yVals, params);
    for (size_t i=0; i<xVals.size(); i++){
        std::cout << expCurve.evaluate(xVals[i]) << " " << yExpVals[i] << std::endl;
        assert(tol >= std::abs(polyCurve.evaluate(xVals[i]) - yExpVals[i])); // should be less than a tolerance
    }


}

/*
Code to read input file. When I have time, I will try to automatically read JSON inputs.
    // Read input file
    std::string file_dir(__FILE__);
    std::string source_dir(file_dir);
    eraseSubStr(source_dir, "main.cpp");
    std::string data_dir(source_dir+"data/");
    std::cout << source_dir << std::endl;
    std::cout << data_dir << std::endl;

    // Short alias for this namespace
    namespace pt = boost::property_tree;
    // Create a root
    pt::ptree KRWIRS;
    // Load the json file in this ptree
    pt::read_json(data_dir+"KRWIRS.json", KRWIRS);
    std::cout << KRWIRS.get("base_date", "None") << std::endl;
    std::cout << KRWIRS.get("tenor_types", "None") << std::endl;
    std::cout << KRWIRS.get("tenors", "None") << std::endl;
    // A vector to allow storing our animals

    std::vector<std::string> str_tenors;
    std::vector<ql::Period> tenors;
    for (pt::ptree::value_type &tenor : KRWIRS.get_child("tenors"))
    {
        // std::cout << tenor.first.data() << std::endl;
        std::cout << tenor.second.data() << std::endl;
        str_tenors.push_back(tenor.second.data());
        tenors.push_back(ql::Period(tenor.second.data()));
    }
*/

#endif /* MAIN_HPP_ */
