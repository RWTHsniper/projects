#include "utils.hpp"


std::vector<std::string> tokenize(std::string s, std::string del)
{
    int start = 0;
    int end = s.find(del);
    std::vector<std::string> output; 
    output.reserve(2);
    while (end != -1) {
        // std::cout << s.substr(start, end - start) << std::endl;
        output.emplace_back(s.substr(start, end - start));
        start = end + del.size();
        end = s.find(del, start);
    }
    output.emplace_back(s.substr(start, end - start));
    return output;
}