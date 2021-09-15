#ifndef SDE_HPP_
#define SDE_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <random>

struct drift{
    size_t lhs_sv;
    double coeff;
    size_t rhs_sv;
    double order;
    drift(std::vector<std::string> p){
        if (p.size() != 4){ std::cout << "Input for drift should be length of 4";};
        lhs_sv = static_cast<size_t>(stoi(p[0]));
        coeff = stod(p[1]);
        rhs_sv = static_cast<size_t>(stoi(p[2]));
        order = stod(p[3]);
    }
    drift(size_t l, double b, size_t r, double o){
        lhs_sv = l;
        coeff = b;
        rhs_sv = r;
        order = o;
    }
};

struct volatility{
    size_t lhs_sv;
    double coeff;
    size_t rhs_sv;
    double order;
    volatility(std::vector<std::string> p){
        if (p.size() != 4){ std::cout << "Input for drift should be length of 4";};
        lhs_sv = static_cast<size_t>(stoi(p[0]));
        coeff = stod(p[1]);
        rhs_sv = static_cast<size_t>(stoi(p[2]));
        order = stod(p[3]);
    }
    volatility(size_t l, double b, size_t r, double o){
        lhs_sv = l;
        coeff = b;
        rhs_sv = r;
        order = o;
    }
};

class SDE {
private:
    double T;
    size_t Nt;
    double dt;
    double sqrt_dt;
    size_t num_paths;
    size_t num_sv; // number of state variables for each model
    std::vector<double> t_vec;
    std::vector<double> x0_vec;
    std::vector<drift> drift_vec;
    std::vector<volatility> volatility_vec;
    std::vector<std::vector<std::vector<double>>> x; // state variables at each step and path (Nt+1, num_sv, num_paths)
    std::vector<std::vector<std::vector<double>>> dW; // Brownian motions at each step and path (Nt, num_sv, num_paths)
    std::vector<std::vector<double>> drift_buffer;
    std::vector<std::vector<double>> volatility_buffer;
public:
	SDE(std::map<std::string, double>& arg_inp_params, std::vector<double>& arg_x0_vec, std::vector<drift>& arg_drift_vec, std::vector<volatility>& arg_volatility_vec);
    void info();
    void compute_drift(size_t ind_t);
    void compute_volatility(size_t ind_t);
    void simulate();
	virtual ~SDE(){}; // empty destructor for a virtual class
};
 

#endif /* SDE_HPP_ */