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

    // get drift info
    // std::vector<size_t, std::tuple<double, size_t, double> > drift_vec;
    std::vector<drift> drift_vec;
    std::vector<volatility> volatility_vec;
    std::vector<double> x0_vec;
    std::vector<correlation> correlation_vec;
    std::map<std::string, double> inp_params;
    for (size_t i=0; i<inp_text.size(); i++){
        std::vector<std::string> &line_str = inp_text[i];
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
    SDE sde = SDE(inp_params, x0_vec, drift_vec, volatility_vec, correlation_vec);
    // sde.compute_drift(0);
    // sde.compute_volatility(0);
    sde.info();
    sde.simulate();
    std::string path  ="../output/";
    sde.write_result(path);

    // drift d_tmp = drift_vec[0];
    // std::cout << d_tmp.lhs_sv << " " << d_tmp.coeff << std::endl;
    // std::map<std::string, double> test;
    // test.insert(std::pair<std::string, double>("test", 2.0));
    // std::cout << test["test"] << std::endl;

    // Cholesky test
    MyTensor test = MyTensor(static_cast<size_t>(3),static_cast<size_t>(3),0.0);
    test.get(0,0)=4.0; test.get(0,1)=12.0; test.get(0,2)=-16.0;
    test.get(1,0)=12.0; test.get(1,1)=37.0; test.get(1,2)=-43.0;
    test.get(2,0)=-16.0; test.get(2,1)=-43.0; test.get(2,2)=98.0;
    test.print(); std::cout << std::endl;
    MyTensor chol_res = test.get_cholesky_lower();
    // *MyTensor<double> chol_res;
    // chol_res.emplace_back(test.get_cholesky_lower());
    
    chol_res.print();

    return 0;
}
