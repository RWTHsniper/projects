#include "stochasticmodel.hpp"

namespace StochasticModel{
    Eigen::VectorXd HaganNF::evolve(const double& t, const Eigen::VectorXd& x, const double& dt, const Eigen::VectorXd& dw) const{
        Eigen::VectorXd xn(x); // state variables at next step (output)
        // Compute drift
        // LGM model does not have any drift
        // Compute volatility
        Eigen::VectorXd vol(nFactor_); 
        for (size_t i=0; i< nFactor_; i++){
            vol[i] = alp[i].evaluate(t);
        }
        for (size_t i=0; i< nFactor_; i++){
            xn[i] += vol[i] * dw[i] * std::sqrt(dt);
        }
        return xn;
    }
}
