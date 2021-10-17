#include "utils.hpp"

std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string remove_comment(const std::string &s, const std::string delimiter){   
    if(s.find(delimiter) == std::string::npos){
        return s;
    }
    else{
        // std::cout << s.find(delimiter) << " # found" << std::endl;
    return s.substr(0, s.find(delimiter)); // token is "scott"
    }
}

std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

std::vector<std::string> tokenize(std::string arg_s, std::string del)
{
    std::string s = rtrim(arg_s);
    int start = 0;
    int end = s.find(del);
    std::vector<std::string> output; 
    output.reserve(2);
    while (end != -1) {
        output.emplace_back(s.substr(start, end - start));
        start = end + del.size();
        end = s.find(del, start);
    }
    output.emplace_back(s.substr(start, end - start));
    return output;
}