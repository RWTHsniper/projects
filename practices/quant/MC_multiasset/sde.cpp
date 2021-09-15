#include "sde.hpp"

SDE::SDE(std::map<std::string, double>& arg_inp_params,  std::vector<drift>& arg_drift_vec, std::vector<volatility>& arg_volatility_vec){
    std::cout << "Initialize SDE" << std::endl;
    T = static_cast<double>(arg_inp_params["T"]);
    Nt = static_cast<size_t>(arg_inp_params["Nt"]);
    dt = T/static_cast<double>(Nt);
    num_paths = static_cast<size_t>(arg_inp_params["num_paths"]);
    drift_vec = arg_drift_vec;
    volatility_vec = arg_volatility_vec;
    num_sv = 0;
    for(std::size_t i = 0; i < drift_vec.size(); ++i) {
        if (drift_vec[i].lhs_sv > num_sv){num_sv = drift_vec[i].lhs_sv;}; 
    }
    num_sv += 1;
    x.resize(Nt+1);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i].resize(num_sv);
        for (std::size_t j = 0; j < num_sv; ++j) {
            x[i][j].resize(num_paths, 0.0);
        }
    }
    //  std::fill(v.begin(), v.end(), c);
    drift_buffer.resize(num_sv);
    volatility_buffer.resize(num_sv); 
    for (std::size_t j = 0; j < num_sv; ++j) {
        drift_buffer[j].resize(num_paths, 0.0);
        volatility_buffer[j].resize(num_paths, 0.0);
    }
    std::cout << "Initialization complete " << std::endl;

}

void SDE::compute_drift(size_t ind_t){
    std::cout << "Test compute drift" << std::endl;
    for (size_t ind_drift=0; ind_drift< drift_vec.size(); ind_drift++){
        drift& elem = drift_vec[ind_drift];
        for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
            for (size_t ind_path=0; ind_path<num_paths; ind_path++){
                drift_buffer[elem.lhs_sv][ind_path]  = elem.coeff*std::pow(x[ind_t][elem.rhs_sv][ind_path],elem.order);
                std::cout << "Drift "<< drift_buffer[elem.lhs_sv][ind_path] << std::endl;
            }
        }
    }
}

void SDE::compute_volatility(size_t ind_t){
    std::cout << "Test compute volatility" << std::endl;
    for (size_t ind_volatility=0; ind_volatility< volatility_vec.size(); ind_volatility++){
        volatility& elem = volatility_vec[ind_volatility];
        for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
            for (size_t ind_path=0; ind_path<num_paths; ind_path++){
                volatility_buffer[elem.lhs_sv][ind_path]  = elem.coeff*std::pow(x[ind_t][elem.rhs_sv][ind_path],elem.order);
                std::cout << "Volatility "<< volatility_buffer[elem.lhs_sv][ind_path] << std::endl;
            }
        }
    }
}