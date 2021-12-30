#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <iostream>
#include <string> // stoi, stod
#include <typeinfo>


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
