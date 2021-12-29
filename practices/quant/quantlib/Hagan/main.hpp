#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <iostream>
#include <string> // stoi, stod
#include <typeinfo>


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <ql/quantlib.hpp>
#include <ql/utilities/dataparsers.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>

// #include <ql/processes/blackscholesprocess.hpp>

#include "utils.hpp"

/*
Code to read input file. When I have time, I will try to automatically read JSON inputs.
    // Read input file
    std::string file_dir(__FILE__);
    std::string source_dir(file_dir);
    eraseSubStr(source_dir, "main.cpp");
    std::string data_dir(source_dir+"data/");
    std::cout << source_dir << std::endl;
    std::cout << data_dir << std::endl;

    // Short alias for this namespace
    namespace pt = boost::property_tree;
    // Create a root
    pt::ptree KRWIRS;
    // Load the json file in this ptree
    pt::read_json(data_dir+"KRWIRS.json", KRWIRS);
    std::cout << KRWIRS.get("base_date", "None") << std::endl;
    std::cout << KRWIRS.get("tenor_types", "None") << std::endl;
    std::cout << KRWIRS.get("tenors", "None") << std::endl;
    // A vector to allow storing our animals

    std::vector<std::string> str_tenors;
    std::vector<ql::Period> tenors;
    for (pt::ptree::value_type &tenor : KRWIRS.get_child("tenors"))
    {
        // std::cout << tenor.first.data() << std::endl;
        std::cout << tenor.second.data() << std::endl;
        str_tenors.push_back(tenor.second.data());
        tenors.push_back(ql::Period(tenor.second.data()));
    }
*/

#endif /* MAIN_HPP_ */
