/*
 * utils.hpp
 *
 *  Created on: Sep 15, 2021
 *      Author: jaeyong
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <vector>
#include <iostream>

// Strings
const std::string WHITESPACE = " \n\r\t\f\v";
std::string ltrim(const std::string &s);
std::string rtrim(const std::string &s);
std::string remove_comment(const std::string &s, const std::string delimiter);
std::string trim(const std::string &s);
std::vector<std::string> tokenize(std::string s, std::string del = " ");

#endif /* UTILS_HPP_ */
