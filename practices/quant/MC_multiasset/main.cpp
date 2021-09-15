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
    std::vector<Models> model_vect;
    for (size_t i=0; i<inp_text.size(); i++){
        for (size_t j=0; j<inp_text[i].size(); j++){
            // std::cout << inp_text[i][j] << " ";
        }
        // std::cout << std::endl;
    }
    for (size_t i=0; i<inp_text.size(); i++){
        std::vector<std::string> &line_str = inp_text[i];
        std::cout << line_str[0] << std::endl;
        std::vector<std::string> params = std::vector<std::string>(line_str.begin() + 2, line_str.end());

        if (line_str[0]=="Model"){
            if (line_str[1]=="GBM"){
                model_vect.emplace_back(GBM(stod(params[0]),stod(params[1])));
            }
        }
    }
    // std::cout << vec[0][1] << std::endl;
 
}
