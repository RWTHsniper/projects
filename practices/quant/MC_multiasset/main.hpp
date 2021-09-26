#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <any>
#include <memory> // for std::unique_ptr

#include "utils.hpp"
#include "sde.hpp"
#include "MyTensor.hpp"

void extract_inputs(std::vector<drift>& drift_vec, std::vector<volatility>& volatility_vec,
    std::vector<double>& x0_vec, std::vector<correlation>& correlation_vec, std::vector<std::unique_ptr<Constraint>>& constraint_vec, 
    std::map<std::string, double>& inp_params, const std::vector<std::vector<std::string>>& inp_text){

    for (size_t i=0; i<inp_text.size(); i++){
        const std::vector<std::string> &line_str = inp_text[i];
        // Check whether the line is empty or not
        size_t tot_line_len = 0;
        for (auto const& e : line_str) {
            tot_line_len += e.size();
        }
        if (tot_line_len == 0){
            continue; // skip empty line
        }
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
        else if (key=="constraint"){
            for (size_t j=0; j<params.size();j++){
                auto p = tokenize(params[j],",");
                // constraint_vec.emplace_back(Contraint_factory(p));
                constraint_vec.push_back(Contraint_factory(p));
            }
        }
        else{
            inp_params.insert(std::pair<std::string, double>(key, stod(line_str[1])));
        }
    }    
}


#endif /* MAIN_HPP_ */
