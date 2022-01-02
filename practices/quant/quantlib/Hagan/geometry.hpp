/*
 * geometry.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include <iostream>
#include <Eigen/Dense> // Dense matrices
// #include <unsupported/Eigen/Polynomials>

namespace Geometry{
    // // Abstract classes
    // class Geometry{
    //     public:
    //         Geometry(){
    //             std::cout << "This is an abstract class. Please make an instance of derived classes" << std::endl;
    //             exit(-1);}
    // };
    // class Curve: Geometry{};

    // Class definitions
    class PolynomialCurve{
        /*
        https://gist.github.com/ksjgh/4d5050d0e9afc5fdb5908734335138d0
        */
        public:
            PolynomialCurve(const Eigen::VectorXd& xVals, const Eigen::VectorXd& yVals, const size_t& order): 
                                                    xVals_(xVals), yVals_(yVals), order_(order){
                                                        std::cout << "x " << xVals_ << std::endl;
                                                        std::cout << "y " << yVals_ << std::endl;
                                                    assert(xVals_.size() == yVals_.size());
                                                    assert(order_ >= 1 && order_ <= xVals_.size() - 1);
                                                    A_.resize(xVals_.size(), order_ + 1);
                                                    coeffs_.resize(order_ + 1);
                                                    computeFitting();}; // copy constructor
        double polyEval(double x);
        private:
            size_t order_;
            Eigen::VectorXd xVals_;
            Eigen::VectorXd yVals_;
            Eigen::VectorXd coeffs_; // polynomial coefficients [0th, 1th, 2th, ...]
            Eigen::MatrixXd A_;
            void computeFitting();
    };


}


#endif /* GEOMETRY_HPP_ */
