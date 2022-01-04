#include "stochasticmodel.hpp"

namespace StochasticModel{
    void HaganNF::evolve(Eigen::MatrixXd& xn, const double& t, const Eigen::MatrixXd& x, const double& dt, const Eigen::MatrixXd& dw) const{
        // Eigen::VectorXd xn(x); // state variables at next step (output)
        // Compute drift
        // LGM model does not have any drift
        // Compute volatility
        // Eigen::MatrixXd vol(xn.size());
        // std::cout << "Start stochastic" << std::endl;
        // std::cout << "rows and cols " << x.rows() << " " << x.cols() << std::endl;
        double vol = 0.0;
        double sqrtDt = std::sqrt(dt);
        for (size_t j=0; j< x.cols(); j++){
            for (size_t i=0; i< x.rows(); i++){
                // vol(i,j) = alp[i].evaluate(t);
                vol = alp[i].evaluate(t);
                xn(i,j) = x(i,j) + vol * dw(i,j) * sqrtDt;
            }
        }
        // for (size_t i=0; i< nFactor_; i++){
        //     vol[i] = alp[i].evaluate(t);
        // }
        // for (size_t i=0; i< nFactor_; i++){
        //     xn[i] += vol[i] * dw[i] * std::sqrt(dt);
        // }
    }
    Eigen::VectorXd HaganNF::evolve(const double& t, const Eigen::VectorXd& x, const double& dt, const Eigen::VectorXd& dw) const{
        Eigen::VectorXd xn(x); // state variables at next step (output)
        double sqrtDt = std::sqrt(dt);
        // Compute drift
        // LGM model does not have any drift
        // Compute volatility
        Eigen::VectorXd vol(nFactor_); 
        for (size_t i=0; i< nFactor_; i++){
            vol[i] = alp[i].evaluate(t);
        }
        for (size_t i=0; i< nFactor_; i++){
            xn[i] += vol[i] * dw[i] * sqrtDt;
        }
        return xn;
    }
    std::shared_ptr<Eigen::MatrixXd> HaganNF::computeInterestRate(const double& t_i, const double& dt, const size_t& numPaths, const size_t& numSteps, std::vector<Eigen::MatrixXd>& x){
        double t = t_i + dt; // skip 0th step
        double T = t_i + numSteps * dt;
        std::shared_ptr<Eigen::MatrixXd> r = std::make_shared<Eigen::MatrixXd>(numPaths, numSteps+1);
        r->setZero();
        // Eigen::MatrixXd r(numPaths, numSteps+1); r.setZero();
        // temporary variables for H, H', and zeta
        Eigen::VectorXd H_t(nFactor_);
        Eigen::VectorXd H_t_deriv(nFactor_);
        Eigen::MatrixXd zeta_t(nFactor_, nFactor_);
        for (size_t step=1; step <= numSteps; step++){ // index for steps
            double fwdRate = yieldCurve_->forwardRate(0.0,6.0,Continuous);
            for (size_t sv=0; sv<nFactor_; sv++){ // index for state variable
                H_t[sv] = H_[sv].evaluate(t);
                H_t_deriv[sv] = H_[sv].evalDeriv(t);
                for (size_t sv2=0; sv2<nFactor_; sv2++){ // index for state variable
                    zeta_t(sv, sv2) = zeta_[sv][sv2].evalInt(t_i, t);
                }
            }
            for (size_t path=0; path < numPaths; path++){ // index for paths
                (*r)(path, step) = fwdRate;
                (*r)(path, step) += H_t_deriv.dot(x[step].col(path) + zeta_t * H_t);
            }
            // // dw = lowerMat_ * dWIndep[i];
            // dw = lowerMat_; dw *= dWIndep[i];
            // haganModel.evolve(x[i+1], t, x[i], dt, dw);
            // if (i==numSteps-1) saveData(source_dir+"output/dw.csv", dw);
            t+= dt;
        }
    return r;
    }
}
