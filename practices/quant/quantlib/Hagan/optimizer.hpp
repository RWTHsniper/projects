/*
 * stochasticmodel.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef OPTIMIZER_HPP_
#define OPTIMIZER_HPP_

#include <iostream>

#include <Eigen/Dense> // Dense matrices
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/NonLinearOptimization>

namespace Optimizer{

   // https://stackoverflow.com/questions/18509228/how-to-use-the-eigen-unsupported-levenberg-marquardt-implementation
   // Generic functor
   // See http://eigen.tuxfamily.org/index.php?title=Functors
   // C++ version of a function pointer that stores meta-data about the function
    template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
    struct Functor
    {

    // Information that tells the caller the numeric type (eg. double) and size (input / output dim)
    typedef _Scalar Scalar;
    enum { // Required by numerical differentiation module
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };

    // Tell the caller the matrix sizes associated with the input, output, and jacobian
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    // Local copy of the number of inputs
    int m_inputs, m_values;

    // Two constructors:
    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    // Get methods for users to determine function input and output dimensions
    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

    };



    // https://en.wikipedia.org/wiki/Test_functions_for_optimization
    // Booth Function
    // Implement f(x,y) = (x + 2*y -7)^2 + (2*x + y - 5)^2
    struct BoothFunctor : Functor<double>
    {
        // Simple constructor
        BoothFunctor(): Functor<double>(2,2) {}

        // Implementation of the objective function
        int operator()(const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {
           double x = z(0);   double y = z(1);
            /*
            * Evaluate the Booth function.
            * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
            * of squared terms. The algorithm takes this into account: do not do it yourself.
            * In other words: objFun = sum(fvec(i)^2)
            */
            fvec(0) = x + 2*y - 7;
            fvec(1) = 2*x + y - 5;
            return 0;
        }
    };


    // https://en.wikipedia.org/wiki/Test_functions_for_optimization
    // Himmelblau's Function
    // Implement f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
    struct HimmelblauFunctor : Functor<double>
    {
        // Simple constructor
        HimmelblauFunctor(): Functor<double>(2,2) {}

        // Implementation of the objective function
        int operator()(const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {
            double x = z(0);   double y = z(1);
            /*
            * Evaluate Himmelblau's function.
            * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
            * of squared terms. The algorithm takes this into account: do not do it yourself.
            * In other words: objFun = sum(fvec(i)^2)
            */
            fvec(0) = x * x + y - 11;
            fvec(1) = x + y * y - 7;
            return 0;
        }
    };

    inline void LMReturnStatus(int ret){
        /*
        LevenbergMarquardt method's return status
        */
        switch(ret){
            case Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation  : std::cout << "Too many function evaluations\n";   break;
            case Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall  : std::cout << "Relative reduction is too small\n";   break;
            case Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall  : std::cout << "Relative error is too small\n";   break;
            case Eigen::LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall  : std::cout << "Relative error and reduction are too small\n";   break;
        }
    }

    void testBoothFun();

}
#endif /* OPTIMIZER_HPP_ */
