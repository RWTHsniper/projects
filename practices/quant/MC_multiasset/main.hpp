#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <any>

#include "utils.hpp"
#include "models.hpp"
#include "sde.hpp"
#include "MyTensor.hpp"

void extract_inputs(std::vector<drift>& drift_vec, std::vector<volatility>& volatility_vec,
    std::vector<double>& x0_vec, std::vector<correlation>& correlation_vec, std::map<std::string, double>& inp_params, 
    const std::vector<std::vector<std::string>>& inp_text){

    for (size_t i=0; i<inp_text.size(); i++){
        const std::vector<std::string> &line_str = inp_text[i];
        auto key = line_str[0];
        std::vector<std::string> params = std::vector<std::string>(line_str.begin() + 1, line_str.end());
        if (key=="drift"){
            for (size_t j=0; j<params.size();j++){
                auto p = tokenize(params[j],",");
                drift_vec.emplace_back(drift(p));
                size_t idx = drift_vec.size()-1;
            }
        }
        else if (key=="volatility"){
            for (size_t j=0; j<params.size();j++){
                auto p = tokenize(params[j],",");
                volatility_vec.emplace_back(volatility(p));
            }
        }
        else if (key=="x0"){
            for (size_t j=0; j<params.size();j++){
                x0_vec.emplace_back(stod(params[j]));
            }
        }
        else if (key=="correlation"){
            for (size_t j=0; j<params.size();j++){
                auto p = tokenize(params[j],",");
                correlation_vec.emplace_back(correlation(stoi(p[0]),stoi(p[1]),stod(p[2])));
            }
        }
        else{
            inp_params.insert(std::pair<std::string, double>(key, stod(line_str[1])));
        }
    }    
}


#endif /* MAIN_HPP_ */
