/*
 * utils.hpp
 *
 *  Created on: Sep 15, 2021
 *      Author: jaeyong
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <vector>
#include <string>
#include <cstddef>
#include <iostream>

#include <boost/property_tree/ptree.hpp>

// Strings
const std::string WHITESPACE = " \n\r\t\f\v";
std::string ltrim(const std::string &s);
std::string rtrim(const std::string &s);
std::string remove_comment(const std::string &s, const std::string delimiter);
std::string trim(const std::string &s);
std::vector<std::string> tokenize(std::string s, std::string del = " ");
void eraseSubStr(std::string & mainStr, const std::string & toErase);

// json helpers
template <typename T>
struct my_id_translator
{
    typedef T internal_type;
    typedef T external_type;

    boost::optional<T> get_value(const T &v) { return  v.substr(1, v.size() - 2) ; }
    boost::optional<T> put_value(const T &v) { return '"' + v +'"'; }
};

#endif /* UTILS_HPP_ */
