#include "sde.hpp"

std::unique_ptr<Constraint> Contraint_factory(std::vector<std::string>& arg_params){
    auto& type = arg_params[0];
    if (type=="min"){
        return std::make_unique<constraint_min>(static_cast<size_t>(stoi(arg_params[1])), stod(arg_params[2]));
    }
    else if (type=="max"){
        return std::make_unique<constraint_max>(static_cast<size_t>(stoi(arg_params[1])), stod(arg_params[2]));
    }
    else{
        std::cout << "Wrong constraint type is given " << std::endl;
        exit(-1);
    }
}

SDE::SDE(std::map<std::string, double>& arg_inp_params, std::vector<double>& arg_x0_vec, std::vector<drift>& arg_drift_vec, std::vector<volatility>& arg_volatility_vec \
, std::vector<correlation>& arg_correlation_vec, std::vector<std::unique_ptr<Constraint>>& arg_constraint_vec){
    std::cout << "Initialize SDE" << std::endl;
    bool use_bound_check = false;
    T = static_cast<double>(arg_inp_params["T"]);
    Nt = static_cast<size_t>(arg_inp_params["Nt"]);
    dt = T/static_cast<double>(Nt);
    sqrt_dt = std::sqrt(dt);
    num_paths = static_cast<size_t>(arg_inp_params["num_paths"]);
    drift_vec = arg_drift_vec;
    volatility_vec = arg_volatility_vec;
    for (size_t i=0; i<arg_constraint_vec.size(); i++){
        constraint_vec.emplace_back(std::move(arg_constraint_vec[i]));
    }
    // Compute num_sv from drift_vec
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
    // Initialization for the state variable x
    x.resize(std::vector<size_t>{Nt+1, num_sv, num_paths}); 
    // Set initial value x0 to x
    for (std::size_t j = 0; j < num_sv; ++j) {
        for (std::size_t k=0; k<num_paths;++k){
            x.get(0,j,k) = x0_vec[j];
        }
    }
    // Initialization for Brownian motion
    std::default_random_engine generator{0}; // Seed 0
    std::normal_distribution<double> distribution(0.0,1.0);
    dW_indep.resize(std::vector<size_t>{Nt, num_sv, num_paths});
    for (std::size_t i = 0; i < Nt; ++i) {
        for (std::size_t j = 0; j < num_sv; ++j) {
            for (std::size_t k=0; k<num_paths;++k){
                dW_indep.get(i,j,k) = distribution(generator);
            }
        }
    }
    // Normalize dW_indep by trnasformation
    for (std::size_t i = 0; i < Nt; ++i) {
        for (std::size_t j = 0; j < num_sv; ++j) {
            double m = dW_indep.mean(i,j,0,i,j,num_paths-1);
            double s = dW_indep.std(i,j,0,i,j,num_paths-1);
            for (std::size_t k=0; k<num_paths;++k){
                dW_indep.get(i,j,k) = (dW_indep.get(i,j,k) - m)/s; // almost exact normal distribution
            }
        }
    }
    //  drift and volatility buffer initialization
    drift_buffer.resize(std::vector<size_t>{num_sv, num_paths});
    volatility_buffer.resize(std::vector<size_t>{num_sv, num_paths});
    // Cholesky matrix
    if ((arg_correlation_vec.size()>0) && (num_sv>1)){
        use_cholesky = true;
        MyTensor correlation_mat = MyTensor(num_sv,num_sv,0.0);
        for (size_t i=0; i<num_sv; i++){correlation_mat.get(i,i) = 1.0;}; // set diagonal terms
        // arg_correlation_vec;
        for (auto & elem : arg_correlation_vec) {
            correlation_mat.get(elem.sv_1,elem.sv_2) = elem.corr;
            correlation_mat.get(elem.sv_2,elem.sv_1) = elem.corr;
        }
        std::cout << "Correlaiton matrix" << std::endl;
        correlation_mat.print();
        cholesky_lower = MyTensor(correlation_mat.get_cholesky_lower()); // Assign cholesky_lower
    }
    else{use_cholesky=false;}
    // Compute dW using the Cholesky lower triangular matrix
    dW.resize(std::vector<size_t>{Nt,num_sv, num_paths});
    for (size_t ind_t=0; ind_t<Nt; ind_t++){
        for (size_t ind_path=0; ind_path<num_paths; ind_path++){
            for (size_t i=0; i<num_sv; i++){ // index for dW
                for (size_t j=0; j<num_sv; j++){ // index for dW_indep
                    if (use_cholesky){
                        dW.get(ind_t,i,ind_path) += cholesky_lower.get(i,j)*dW_indep.get(ind_t,j,ind_path);
                    }
                    else{
                        dW.get(ind_t,i,ind_path) += dW_indep.get(ind_t,j,ind_path);
                    }
                }
            }
        }
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
    cholesky_lower.print(); std::cout << std::endl; 
}

void SDE::compute_drift(size_t ind_t){
    // std::cout << "Test compute drift" << std::endl;
    drift_buffer = static_cast<double>(0.0); // initialization
    for (size_t ind_drift=0; ind_drift< drift_vec.size(); ind_drift++){
        drift& elem = drift_vec[ind_drift];
        for (size_t ind_path=0; ind_path<num_paths; ind_path++){
            drift_buffer.get(elem.lhs_sv,ind_path) += elem.compute(x.get(ind_t, elem.rhs_sv, ind_path),dt);
        }
    }
}

void SDE::compute_volatility(size_t ind_t){
    // std::cout << "Test compute volatility" << std::endl;
    volatility_buffer = double(0.0); // initialization
    for (size_t ind_volatility=0; ind_volatility< volatility_vec.size(); ind_volatility++){
        volatility& elem = volatility_vec[ind_volatility];
        for (size_t ind_path=0; ind_path<num_paths; ind_path++){
            // Use Milstein scheme
            double dW_t = dW.get(ind_t,elem.lhs_sv,ind_path)*sqrt_dt;
            volatility_buffer.get(elem.lhs_sv,ind_path) += elem.compute_milstein(x.get(ind_t, elem.rhs_sv, ind_path),dW_t,dt); // conventional
        }
    }
}

void SDE::simulate(){
    // omp_set_num_threads(num_threads);
    // #pragma omp parallel{
    //     int tid = omp_get_thread_num();
    // std::cout << "I am " << tid << std::endl;
    // }


    std::cout << "Start simulation" << std::endl;
    for (size_t ind_t=0; ind_t<t_vec.size()-1;ind_t++){
        this->compute_drift(ind_t); // compute drift
        this->compute_volatility(ind_t); // compute volatility
        // accumulate
        for (size_t ind_sv=0; ind_sv<num_sv; ind_sv++){
            for (size_t ind_path=0; ind_path<num_paths; ind_path++){
                x.get(ind_t+1,ind_sv,ind_path) = x.get(ind_t,ind_sv,ind_path);
                x.get(ind_t+1,ind_sv,ind_path) += drift_buffer.get(ind_sv,ind_path);
                x.get(ind_t+1,ind_sv,ind_path) += volatility_buffer.get(ind_sv,ind_path);
            }
        }
        if (constraint_vec.size() > 0){ // Apply constraints when at least one constraint object is declared
            for (size_t i_c=0; i_c< constraint_vec.size(); i_c++){ // loop over and apply each constraint object
                const size_t& ind_sv = constraint_vec[i_c]->get_index(); // Index to be constrained
                for (size_t ind_path=0; ind_path<num_paths; ind_path++){
                    x.get(ind_t+1,ind_sv,ind_path) = constraint_vec[i_c]->compute(x.get(ind_t+1,ind_sv,ind_path));
                }
            }
        }
    }
    std::cout << "Simulation is complete" << std::endl;
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
                myFile << x.get(ind_t,ind_sv,ind_path); 
                // myFile << x[ind_t][ind_sv][ind_path]; 
                if (ind_t<t_vec.size()-1){myFile << ",";};
            }
            myFile << "\n";
        }
        // Close the file
        myFile.close();
    }

}

void SDE::print_result(){
    std::cout << "Print statistics at the maturity" << std::endl;
    for (size_t j=0; j<num_sv; j++){
        std::cout << j << "-th state variable" << std::endl;
        std::cout << "mean: " << x.mean(Nt, j, 0, Nt, j, num_paths-1) << " ";
        std::cout << "std: " << x.std(Nt, j, 0, Nt, j, num_paths-1);
        std::cout << std::endl;
    }
}