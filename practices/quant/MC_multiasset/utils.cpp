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