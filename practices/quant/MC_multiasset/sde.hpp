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
#include <omp.h> // openmp -fopenmp
#include <chrono> // measure time

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
    virtual void compute_omp(MyTensor<double>& x, const size_t& ind_t, const size_t& num_paths){
        std::cout << "parent compute_omp" << std::endl;
    };
};

class constraint_max: public Constraint {
private:
	double max_value; // minimum value for a constraint
public:
    constraint_max(const size_t arg_ind, double arg_value){
        ind = arg_ind;
        max_value = arg_value;
    };
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
	double min_value; // minimum value for a constraint
public:
    constraint_min(const size_t arg_ind, double arg_value){
        ind = arg_ind;
        min_value = arg_value;
    };
    const size_t& get_index(){return ind;};
    bool check_index(const size_t& arg_ind){
        return (ind==arg_ind);
    }
    double compute(const double& x) override{
        // std::cout << "x "<<x<< " val "<<min_value << " returning " << std::max(x,min_value) << std::endl;
        return std::max(x, min_value); // returned value should be larger than min_value
    };
    void compute_omp(MyTensor<double>& x, const size_t& ind_t, const size_t& num_paths) override{ // slower than single thread
        // std::cout << "Start compute_omp" << std::endl;
        // auto time_start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for schedule(static)
        for (size_t ind_path=0; ind_path<num_paths; ind_path++){
            // x.get(ind_t+1,ind,ind_path) = compute(x.get(ind_t+1,ind,ind_path));
            x.get(ind_t+1,ind,ind_path) = std::max(x.get(ind_t+1,ind,ind_path), min_value);
        }
    //     auto time_end = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);
    //     std::cout << "Time for running compute_omp " << (duration.count()) << " micro sec" << std::endl;
    };
	~constraint_min(){};
};

std::unique_ptr<Constraint> Contraint_factory(std::vector<std::string>& arg_params);

class PoissonProcess {
    private:
        double coeff;
        double intensity; // intensity
        double mean; // intensity*dt
        double dt;
    public:
        size_t lhs_sv;
        PoissonProcess(std::vector<std::string>& p){
            if (p.size() != 3){ std::cout << "Input for drift should be length of 3";};
            lhs_sv = static_cast<size_t>(stoi(p[0]));
            coeff = stod(p[1]);
            intensity = stod(p[2]);
        }
        PoissonProcess(const size_t l, const double& c, const double& in) : lhs_sv(l), coeff(c), intensity(in) {};
        void set_dt(const double& arg_dt){
            dt = arg_dt;
        }
        double get_intensity(){return intensity;}
        double compute(){
            return coeff;
        }
        double compute(const double& arg_dt){
            mean = intensity*arg_dt;
            return coeff;
        }
        ~PoissonProcess(){};
};

class SDE {
private:
    double T;
    size_t Nt;
    double dt;
    double sqrt_dt;
    size_t num_paths;
    size_t num_sv; // number of state variables for each model
    size_t num_poisson;
    size_t num_threads; // openmp threads
    std::vector<double> t_vec;
    std::vector<double> x0_vec;
    std::vector<drift> drift_vec;
    std::vector<volatility> volatility_vec;
    std::vector<PoissonProcess> poisson_vec;
    std::vector<std::unique_ptr<Constraint>> constraint_vec;
    MyTensor<double> cholesky_lower; // Cholesky lower matrix
    bool use_cholesky;
    MyTensor<double> x; // state variables at each step and path (Nt+1, num_sv, num_paths)
    MyTensor<double> dW; // Brownian motions at each step and path (Nt, num_sv, num_paths)
    MyTensor<double> dW_indep; // Independent Brownian motions at each step and path (Nt, num_sv, num_paths)
    MyTensor<double> dN; // Poisson process at each step and path (Nt, num_sv, num_paths)
    MyTensor<double> drift_buffer;
    MyTensor<double> volatility_buffer;
public:
	SDE(std::map<std::string, double>& arg_inp_params, std::vector<double>& arg_x0_vec, std::vector<drift>& arg_drift_vec, std::vector<volatility>& arg_volatility_vec \
    , std::vector<correlation>& arg_correlation_vec, std::vector<PoissonProcess>& arg_poisson_vec, std::vector<std::unique_ptr<Constraint>>& arg_constraint_vec);
    void info();
    void compute_drift(size_t ind_t);
    void compute_volatility(size_t ind_t);
    void simulate();
    void write_result(std::string& filename);
    void print_result();
	virtual ~SDE(){}; // empty destructor for a virtual class
};
 

#endif /* SDE_HPP_ */