#include "optimizer.hpp"

void Optimizer::testBoothFun() {
    std::cout << "Testing the Booth function..." << std::endl;
    Eigen::VectorXd zInit(2); zInit << 1.87, 2.032;
    std::cout << "zInit: " << zInit.transpose() << std::endl;
    Eigen::VectorXd zSoln(2); zSoln << 1.0, 3.0;
    std::cout << "zSoln: " << zSoln.transpose() << std::endl;

    BoothFunctor functor;
    Eigen::NumericalDiff<BoothFunctor> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BoothFunctor>,double> lm(numDiff);
    lm.parameters.maxfev = 1000;
    lm.parameters.xtol = 1.0e-10;
    std::cout << "max fun eval: " << lm.parameters.maxfev << std::endl;
    std::cout << "x tol: " << lm.parameters.xtol << std::endl;

    Eigen::VectorXd z = zInit;
    int ret = lm.minimize(z);
    std::cout << "iter count: " << lm.iter << std::endl;
    std::cout << "return status: " << ret << std::endl;
    std::cout << "zSolver: " << z.transpose() << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}