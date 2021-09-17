#ifndef SDE_HPP_
#define SDE_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <random>
#include <memory> // for std::unique_ptr
#include <utility> // for std::move

#include "MyTensor.hpp"

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
    double compute(double& sv){
        return coeff*std::pow(sv,order);
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
    double compute(double& sv){
        return coeff*std::pow(sv,order);
    }
};

struct correlation{
    size_t sv_1;
    size_t sv_2;
    double corr;
    correlation(size_t asv_1, size_t asv_2, double acorr):sv_1(asv_1),sv_2(asv_2),corr(acorr) {};
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
    std::unique_ptr<MyTensor<double>> cholesky_ptr;
    bool use_cholesky;
    std::vector<std::vector<std::vector<double>>> x; // state variables at each step and path (Nt+1, num_sv, num_paths)
    std::unique_ptr<MyTensor<double>> dW; // Brownian motions at each step and path (Nt, num_sv, num_paths)
    std::unique_ptr<MyTensor<double>> dW_indep; // Independent Brownian motions at each step and path (Nt, num_sv, num_paths)
    // std::vector<std::vector<std::vector<double>>> dW_indep; // Independent Brownian motions at each step and path (Nt, num_sv, num_paths)
    std::vector<std::vector<double>> drift_buffer;
    std::vector<std::vector<double>> volatility_buffer;
public:
	SDE(std::map<std::string, double>& arg_inp_params, std::vector<double>& arg_x0_vec, std::vector<drift>& arg_drift_vec, std::vector<volatility>& arg_volatility_vec \
    , std::vector<correlation>& arg_correlation_vec);
    void info();
    void compute_drift(size_t ind_t);
    void compute_volatility(size_t ind_t);
    void simulate();
    void write_result(std::string& filename);
	virtual ~SDE(){}; // empty destructor for a virtual class
};
 

#endif /* SDE_HPP_ */