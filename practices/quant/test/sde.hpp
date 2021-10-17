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
#include <omp.h> // openmp

#include "MyTensor.hpp"

struct drift{
    size_t lhs_sv;
    double coeff;
    size_t rhs_sv;
    double order;
    drift(std::vector<std::string>& p){
        if (p.size() != 4){ std::cout << "Input for drift should be length of 4";};
        lhs_sv = static_cast<size_t>(stoi(p[0]));
        coeff = stod(p[1]);
        rhs_sv = static_cast<size_t>(stoi(p[2]));
        order = stod(p[3]);
    }
    drift(const size_t l, const double b, const size_t r, const double o){
        lhs_sv = l;
        coeff = b;
        rhs_sv = r;
        order = o;
    }
    double compute(double& sv, double& dt){
        return coeff*std::pow(sv,order)*dt;
    }
};

struct volatility{
    size_t lhs_sv;
    double coeff;
    size_t rhs_sv;
    double order;
    volatility(std::vector<std::string>& p){
        if (p.size() != 4){ std::cout << "Input for drift should be length of 4";};
        lhs_sv = static_cast<size_t>(stoi(p[0]));
        coeff = stod(p[1]);
        rhs_sv = static_cast<size_t>(stoi(p[2]));
        order = stod(p[3]);
    }
    volatility(const size_t l, const double b, const size_t r, const double o){
        lhs_sv = l;
        coeff = b;
        rhs_sv = r;
        order = o;
    }
    double compute(const double& sv, const double& dt){
        return coeff*std::pow(sv,order)*dt;
    }
    double compute_milstein(const double& sv, const double& dW_t, const double& dt){
        double ret;
        if(order==0){
            ret = coeff*dW_t;
        }
        else{
            ret = coeff*std::pow(sv,order)*dW_t + 0.5*std::pow(coeff,2)*order*std::pow(sv,2*order-1)*(std::pow(dW_t,2)-dt);
        }
        return ret;
    }    
};

struct correlation{
    size_t sv_1;
    size_t sv_2;
    double corr;
    correlation(size_t asv_1, size_t asv_2, double acorr):sv_1(asv_1),sv_2(asv_2),corr(acorr) {};
};

class Constraint {
protected:
    size_t ind;
public:
    virtual bool check_index(const size_t& arg_ind){return true;};
    virtual const size_t& get_index(){return ind;};
    virtual double compute(const double& x){
        std::cout << "parent compute" << std::endl;
        return 0.0;};
};

class constraint_max: public Constraint {
private:
    size_t ind; // Constraint on ind-th state variable
	double max_value; // minimum value for a constraint
public:
    constraint_max(const size_t arg_ind, double arg_value):ind(arg_ind),max_value(arg_value){};
    const size_t& get_index(){return ind;};
    bool check_index(const size_t& arg_ind){
        return (ind==arg_ind);
    }
    double compute(const double& x){
        return std::min(x, max_value); // returned value should be less than min_value
    };
	~constraint_max(){};
};

class constraint_min: public Constraint {
private:
    size_t ind; // Constraint on ind-th state variable
	double min_value; // minimum value for a constraint
public:
    constraint_min(const size_t arg_ind, double arg_value):ind(arg_ind),min_value(arg_value){};
    const size_t& get_index(){return ind;};
    bool check_index(const size_t& arg_ind){
        return (ind==arg_ind);
    }
    double compute(const double& x) {
        // std::cout << "x "<<x<< " val "<<min_value << " returning " << std::max(x,min_value) << std::endl;
        return std::max(x, min_value); // returned value should be larger than min_value
    };
	~constraint_min(){};
};

std::unique_ptr<Constraint> Contraint_factory(std::vector<std::string>& arg_params);

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
    std::vector<std::unique_ptr<Constraint>> constraint_vec;
    MyTensor<double> cholesky_lower; // Cholesky lower matrix
    bool use_cholesky;
    MyTensor<double> x; // state variables at each step and path (Nt+1, num_sv, num_paths)
    MyTensor<double> dW; // Brownian motions at each step and path (Nt, num_sv, num_paths)
    MyTensor<double> dW_indep; // Independent Brownian motions at each step and path (Nt, num_sv, num_paths)
    MyTensor<double> drift_buffer;
    MyTensor<double> volatility_buffer;
public:
	SDE(std::map<std::string, double>& arg_inp_params, std::vector<double>& arg_x0_vec, std::vector<drift>& arg_drift_vec, std::vector<volatility>& arg_volatility_vec \
    , std::vector<correlation>& arg_correlation_vec, std::vector<std::unique_ptr<Constraint>>& arg_constraint_vec);
    void info();
    void compute_drift(size_t ind_t);
    void compute_volatility(size_t ind_t);
    void simulate();
    void write_result(std::string& filename);
    void print_result();
	virtual ~SDE(){}; // empty destructor for a virtual class
};
 

#endif /* SDE_HPP_ */