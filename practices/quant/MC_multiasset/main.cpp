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
        line = remove_comment(line, "#"); // skip comments
        inp_text.emplace_back(tokenize(line));
    }
    // close the file stream
    inFile.close();
    // get drift info
    std::vector<drift> drift_vec;
    std::vector<volatility> volatility_vec;
    std::vector<double> x0_vec;
    std::vector<correlation> correlation_vec;
    std::vector<std::unique_ptr<Constraint>> constraint_vec;
    std::map<std::string, double> inp_params;

    extract_inputs(drift_vec, volatility_vec, x0_vec, correlation_vec, constraint_vec, inp_params, inp_text);
    SDE sde = SDE(inp_params, x0_vec, drift_vec, volatility_vec, correlation_vec, constraint_vec);
    // sde.compute_drift(0);
    // sde.compute_volatility(0);
    sde.info();
    sde.simulate();
    std::string path  ="../output/";
    sde.write_result(path);
    sde.print_result();

    if(false){ // Playground for testing classes
        MyTensor<double> test = MyTensor<double>(3,3,3,2.0,true);
        // std::cout << test.mean(0,1,1,1,1,2) << std::endl;
        test.get(1,1,1) = 10.0; test.get(1,1,2) = 30.0;
        std::cout << test.get(1,1,2) << std::endl;
        std::cout << test.get(std::vector<std::size_t>{1,1,2}) << std::endl;
        test.get(std::vector<std::size_t>{2,1,2}) = 25;
        // std::cout << test.get(std::vector<std::size_t>{2,3,2}) << std::endl;

        std::cout << test.mean(1,1,1,1,1,2) << std::endl;
        std::cout << test.mean(std::vector<std::size_t>{1,1,1},std::vector<std::size_t>{1,1,2}) << std::endl;
        std::cout << test.std(1,1,1,1,1,2) << std::endl;
        std::cout << test.std(std::vector<std::size_t>{1,1,1},std::vector<std::size_t>{1,1,2}) << std::endl;
    }

    return 0;
}
