#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <iostream>
#include <fstream>
#include <string> // stoi, stod
#include <typeinfo>
#include <utility> // move
#include <memory> // smart pointer

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
#include "simulation.hpp"
#include "optimizer.hpp"

using namespace QuantLib;


void testModules(){
    Model::testModule();
    std::cout << "All of the tests for modules are complete!" << std::endl;
    exit(-1); // terminate program
}

void writeOutputJSON(const std::string& filePath, std::shared_ptr<std::vector<Period>> swaptionExpiry, std::shared_ptr<std::vector<Period>> swaptionTenor,
                    Eigen::MatrixXd computedIvol){

    // Write output as JSON
    std::string jsonMsg;
    jsonMsg = "{\"expiry_size\": " + std::to_string(swaptionExpiry->size()) + ",\n";
    jsonMsg += "\"tenor_size\": " + std::to_string(swaptionTenor->size()) + ",\n";
    jsonMsg += "\"expiry\": ["; 
    for (size_t i=0; i < swaptionExpiry->size(); i++){
        jsonMsg += std::to_string(ql::years((*swaptionExpiry)[i]));
        if (i == swaptionExpiry->size()-1) jsonMsg += "],\n";
        else jsonMsg += ", \n";
    }
    jsonMsg += "\"tenor\": ["; 
    for (size_t i=0; i < swaptionTenor->size(); i++){
        jsonMsg += std::to_string(ql::years((*swaptionTenor)[i]));
        if (i == swaptionTenor->size()-1) jsonMsg += "],\n";
        else jsonMsg += ", \n";
    }
    jsonMsg += "\"quote\": ["; 
    // row major
    for (size_t i=0; i<swaptionExpiry->size(); i++){
        for (size_t j=0; j<swaptionTenor->size(); j++){
            jsonMsg += std::to_string(computedIvol(i,j));
            if ((i == swaptionExpiry->size()-1) &&(j == swaptionTenor->size()-1)) jsonMsg += "]\n";
            else jsonMsg += ", \n";
        }
    }
    jsonMsg += "}";
    std::ofstream myfile;
    myfile.open (filePath);
    myfile << jsonMsg;
    myfile.close();    


}

#endif /* MAIN_HPP_ */
