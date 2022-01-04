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
}
