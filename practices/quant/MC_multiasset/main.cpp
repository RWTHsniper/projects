// #include <iostream>

#include "main.hpp"

int main(int, char**) {
    std::ifstream inFile;
    // open the file stream
    inFile.open("../input/test.txt"); 
    // check if opening a file failed
    if (inFile.fail()) {
        std::cerr << "Error opeing a file" << std::endl;
        inFile.close();
        exit(1);
    }
    std::string line;
    std::vector<std::vector<std::string> > inp_text;
    while (std::getline(inFile, line))
    {
        inp_text.emplace_back(tokenize(line));
    }
    // close the file stream
    inFile.close();

    // Print contents in input file
    std::vector<Models> model_vec;
    for (size_t i=0; i<inp_text.size(); i++){
        for (size_t j=0; j<inp_text[i].size(); j++){
            std::cout << inp_text[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // get drift info
    // std::vector<size_t, std::tuple<double, size_t, double> > drift_vec;
    std::vector<drift> drift_vec;
    std::vector<volatility> volatility_vec;
    std::map<std::string, double> inp_params;
    for (size_t i=0; i<inp_text.size(); i++){
        std::vector<std::string> &line_str = inp_text[i];
        auto key = line_str[0];
        std::cout << key << std::endl;
        if (key=="drift"){
            std::vector<std::string> params = std::vector<std::string>(line_str.begin() + 1, line_str.end());
            for (size_t j=0; j<params.size();j++){
                auto p = tokenize(params[j],",");
                drift_vec.emplace_back(drift(p));
                size_t idx = drift_vec.size()-1;
                std::cout<< drift_vec[idx].lhs_sv << " ";
                std::cout<< drift_vec[idx].coeff << " ";
                std::cout<< drift_vec[idx].rhs_sv << " ";
                std::cout<< drift_vec[idx].order << std::endl;
            }
        }
        else if (line_str[0]=="volatility"){
            std::vector<std::string> params = std::vector<std::string>(line_str.begin() + 1, line_str.end());
            for (size_t j=0; j<params.size();j++){
                auto p = tokenize(params[j],",");
                volatility_vec.emplace_back(volatility(p));
                size_t idx = volatility_vec.size()-1;
                std::cout<< volatility_vec[idx].lhs_sv << " ";
                std::cout<< volatility_vec[idx].coeff << " ";
                std::cout<< volatility_vec[idx].rhs_sv << " ";
                std::cout<< volatility_vec[idx].order << std::endl;
            }
        }
        else{
            inp_params.insert(std::pair<std::string, double>(key, stod(line_str[1])));
        }

    }
    std::cout << "Input parameters for simulation" << std::endl;
    for (auto const &pair: inp_params) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
    SDE sde = SDE(inp_params, drift_vec, volatility_vec);
    sde.compute_drift(0);
    sde.compute_volatility(0);
    // drift d_tmp = drift_vec[0];
    // std::cout << d_tmp.lhs_sv << " " << d_tmp.coeff << std::endl;

    /*
    for (size_t i=0; i<inp_text.size(); i++){
        std::vector<std::string> &line_str = inp_text[i];
        std::cout << line_str[0] << std::endl;
        std::vector<std::string> params = std::vector<std::string>(line_str.begin() + 2, line_str.end());

        if (line_str[0]=="Model"){
            if (line_str[1]=="GBM"){
                model_vec.emplace_back(GBM(stod(params[0]),stod(params[1])));
            }
        }
    }
    */
    // std::cout << vec[0][1] << std::endl;
 
}
