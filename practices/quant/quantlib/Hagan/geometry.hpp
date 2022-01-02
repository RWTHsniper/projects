/*
 * geometry.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include <Eigen/Dense> // Dense matrices
// #include <unsupported/Eigen/Polynomials>

// Abstract classes
class Geometry{};
class Curve: Geometry{};

// Class definitions
class PolynomialCurve{
    public:
        PolynomialCurve(const Eigen::VectorXd& xVals, const Eigen::VectorXd& yVals): xVals_(xVals), yVals_(yVals){}; // copy constructor
    private:
        Eigen::VectorXd xVals_;
        Eigen::VectorXd yVals_;
};

#endif /* GEOMETRY_HPP_ */
