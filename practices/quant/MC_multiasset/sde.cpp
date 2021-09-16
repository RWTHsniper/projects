#include "sde.hpp"

SDE::SDE(std::map<std::string, double>& arg_inp_params, std::vector<double>& arg_x0_vec, std::vector<drift>& arg_drift_vec, std::vector<volatility>& arg_volatility_vec \
, std::vector<correlation>& arg_correlation_vec){
    std::cout << "Initialize SDE" << std::endl;
    T = static_cast<double>(arg_inp_params["T"]);
    Nt = static_cast<size_t>(arg_inp_params["Nt"]);
    dt = T/static_cast<double>(Nt);
    sqrt_dt = std::sqrt(dt);
    num_paths = static_cast<size_t>(arg_inp_params["num_paths"]);
    drift_vec = arg_drift_vec;
    volatility_vec = arg_volatility_vec;
    num_sv = 0;
    for(std::size_t i = 0; i < drift_vec.size(); ++i) {
        if (drift_vec[i].lhs_sv > num_sv){num_sv = drift_vec[i].lhs_sv;}; 
    }
    num_sv += 1;
    x0_vec = arg_x0_vec;
    if (num_sv != x0_vec.size()){
        std::cout << "Mismatch between num_sv and x0_vec.size()!" << std::endl;
    }
    t_vec.resize(Nt+1,0.0);
    for (size_t i=0; i<t_vec.size();i++){t_vec[i] = i*dt;};
    x.resize(Nt+1);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i].resize(num_sv);
        for (std::size_t j = 0; j < num_sv; ++j) {
            x[i][j].resize(num_paths, 0.0);
        }
    }
    dW.resize(Nt);
    for (std::size_t i = 0; i < dW.size(); ++i) {
        dW[i].resize(num_sv);
        for (std::size_t j = 0; j < num_sv; ++j) {
            dW[i][j].resize(num_paths);
            for (std::size_t k=0; k<num_paths;++k){
                dW[i][j][k] = distribution(generator);
            }
        }
    }
    // Set initial value x0 to x
    for (std::size_t j = 0; j < num_sv; ++j) {
        for (std::size_t k=0; k<num_paths;++k){
            x[0][j][k] = x0_vec[j];
        }
    }
    //  std::fill(v.begin(), v.end(), c);
    drift_buffer.resize(num_sv);
    volatility_buffer.resize(num_sv); 
    for (std::size_t j = 0; j < num_sv; ++j) {
        drift_buffer[j].resize(num_paths, 0.0);
        volatility_buffer[j].resize(num_paths, 0.0);
    }

    // Cholesky matrix
    if ((arg_correlation_vec.size()>0) && (num_sv>1)){
        MyTensor correlation_mat = MyTensor(num_sv,num_sv,0.0);
        for (size_t i=0; i<num_sv; i++){correlation_mat.get(i,i) = 1.0;}; // set diagonal terms
        // arg_correlation_vec;
        for (auto & elem : arg_correlation_vec) {
            correlation_mat.get(elem.sv_1,elem.sv_2) = elem.corr;
            correlation_mat.get(elem.sv_2,elem.sv_1) = elem.corr;
        }
        std::cout << "Correlaiton matrix" << std::endl;
        correlation_mat.print();
        cholesky_mat.emplace_back(correlation_mat.get_cholesky_lower());
        cholesky_mat[0].print();
    }


    std::cout << "Initialization complete " << std::endl;


}

void SDE::info(){
    std::cout << "Simulation Information" << std::endl;
    std::cout << "Maturity " << T << std::endl;
    std::cout << "Number of time steps " << Nt << std::endl;
    std::cout << "Time step size " << dt << std::endl; 
    std::cout << "Number of paths " << num_paths << std::endl; 
    std::cout << "x0: " ;  for (double x0: x0_vec){std::cout<<x0<<" ";}; std::cout << std::endl; 
    std::cout << "Cholesky lower matrix" << std::endl;
    cholesky_mat[0].print();
    
    std::cout << std::endl; 
}

void SDE::compute_drift(size_t ind_t){
    // std::cout << "Test compute drift" << std::endl;
    for (size_t ind_drift=0; ind_drift< drift_vec.size(); ind_drift++){
        drift& elem = drift_vec[ind_drift];
        for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
            for (size_t ind_path=0; ind_path<num_paths; ind_path++){
                drift_buffer[elem.lhs_sv][ind_path]  = elem.compute(x[ind_t][elem.rhs_sv][ind_path]);
                // assert(elem.compute(x[ind_t][elem.rhs_sv][ind_path]) == elem.coeff*std::pow(x[ind_t][elem.rhs_sv][ind_path],elem.order));
                // drift_buffer[elem.lhs_sv][ind_path]  = elem.coeff*std::pow(x[ind_t][elem.rhs_sv][ind_path],elem.order);
                // std::cout << "Drift "<< drift_buffer[elem.lhs_sv][ind_path] << std::endl;
            }
        }
    }
}

void SDE::compute_volatility(size_t ind_t){
    // std::cout << "Test compute volatility" << std::endl;
    for (size_t ind_volatility=0; ind_volatility< volatility_vec.size(); ind_volatility++){
        volatility& elem = volatility_vec[ind_volatility];
        for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
            for (size_t ind_path=0; ind_path<num_paths; ind_path++){
                volatility_buffer[elem.lhs_sv][ind_path]  = elem.compute(x[ind_t][elem.rhs_sv][ind_path]);
                // assert(elem.compute(x[ind_t][elem.rhs_sv][ind_path]) == elem.coeff*std::pow(x[ind_t][elem.rhs_sv][ind_path],elem.order));
                // volatility_buffer[elem.lhs_sv][ind_path]  = elem.coeff*std::pow(x[ind_t][elem.rhs_sv][ind_path],elem.order);
                // std::cout << "Volatility "<< volatility_buffer[elem.lhs_sv][ind_path] << std::endl;
            }
        }
    }
}

void SDE::simulate(){
    std::cout << "Start simulation" << std::endl;
    for (size_t ind_t=0; ind_t<t_vec.size()-1;ind_t++){
        this->compute_drift(ind_t); // compute drift
        this->compute_volatility(ind_t); // compute volatility
        // accumulate
        for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
            for (size_t ind_path=0; ind_path<num_paths; ind_path++){
                x[ind_t+1][ind_sv][ind_path] = x[ind_t][ind_sv][ind_path];
                x[ind_t+1][ind_sv][ind_path] += drift_buffer[ind_sv][ind_path]*dt;
                x[ind_t+1][ind_sv][ind_path] += volatility_buffer[ind_sv][ind_path]*dW[ind_t][ind_sv][ind_path]*sqrt_dt;
            }
        }
    }
    // Print result
        // for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
        //     std::cout << ind_sv << "-th state variable's paths" << std::endl;
        //     for (size_t ind_path=0; ind_path<num_paths; ind_path++){
        //         std::cout << x[Nt][ind_sv][ind_path] << std::endl;
        //     }
        // }
}

void SDE::write_result(std::string& path){

    // Write output of x
    for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
        std::string colname = "x"+std::to_string(ind_sv);
        std::string filename = path+colname+".csv";
        // Create an output filestream object
        std::ofstream myFile(filename);
        
        // Send data to the stream
        for (size_t ind_path=0; ind_path<num_paths; ind_path++){
            for (size_t ind_t=0; ind_t<t_vec.size(); ind_t++){
                myFile << x[ind_t][ind_sv][ind_path]; 
                if (ind_t<t_vec.size()-1){myFile << ",";};
            }
            myFile << "\n";
        }
        // Close the file
        myFile.close();
    }

}