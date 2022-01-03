#include "stochasticmodel.hpp"

namespace StochasticModel{
    Eigen::VectorXd HaganNF::evolve(double t, Eigen::VectorXd x, double dt, Eigen::VectorXd dw) const{
        Eigen::VectorXd xn(nFactor_); // state variables at next step (output)
        return xn;
    }
}
