/*
 * stochasticmodel.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef STOCHASTICMODEL_HPP_
#define STOCHASTICMODEL_HPP_

#include <iostream>

#include <ql/quantlib.hpp>
#include <ql/time/calendar.hpp>
#include <ql/pricingengine.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>
#include <Eigen/Dense> // Dense matrices

#include "model.hpp"

using namespace QuantLib;

namespace StochasticModel{
   class HaganNF{
      public:
         HaganNF(ext::shared_ptr<YieldTermStructure>& yieldCurve, const size_t& nFactor, Eigen::MatrixXd& corrMat): yieldCurve_(yieldCurve),
                                 nFactor_(nFactor), corrMat_(corrMat){
                                    H_.resize(nFactor);
                                    zeta_.resize(nFactor, nFactor);
                                    if (!corrMat_.isApprox(corrMat_.transpose())){
                                       throw std::runtime_error("Input correlation matrix is not symmetric!");
                                    }
                                    // Apply Cholesky decomposition for the correlation matrix
                                    lowerMat_ = Eigen::MatrixXd(corrMat_.llt().matrixL());
                                    if (!corrMat_.isApprox(lowerMat_*lowerMat_.transpose())){
                                       std::cout << "Correlation matrix " << corrMat_ << std::endl;
                                       std::cout << "Lower matrix " << lowerMat_ << std::endl;
                                       throw std::runtime_error("Cholesky decomposition for the correlation matrix is failed!");
                                    }
                                    Eigen::VectorXd coeffs(2); coeffs << 0.0,1.0;
                                    size_t order = 1;
                                    alp.reserve(nFactor_);
                                    for (size_t i=0; i<nFactor_; i++){alp.emplace_back(order, coeffs);};
         }
         Eigen::MatrixXd getLowerMat() const {return lowerMat_;}
         Eigen::VectorXd evolve(double t, Eigen::VectorXd x, double dt, Eigen::VectorXd dw) const;

      private:
         size_t nFactor_;
         std::vector<Model::PolyFunc> alp; // alpha
         Eigen::VectorXd H_; // How to implement H_ in the framework?
         Eigen::MatrixXd zeta_;
         Eigen::MatrixXd corrMat_; // correlation matrices b.t.w. factors
         Eigen::MatrixXd lowerMat_; // Lower part of the Cholesky decomposition of corrMat_
         ext::shared_ptr<YieldTermStructure> yieldCurve_; // computes discount factor and forward rate


   };
}


/*
// Base class
using namespace std;
class Shape {
   public:
      void setWidth(int w) {
         width = w;
      }
      void setHeight(int h) {
         height = h;
      }
      template <typename F> void test(F arg0){
        std::cout << "test function in the Shape class must be implmented in children classes!" << std::endl;
        exit(-1);
      }
      
   protected:
      int width;
      int height;
};

    //   struct drift{
    //     size_t lhs_sv;
    //     double coeff;
    //     size_t rhs_sv;
    //     double order;
    //   };    

// Derived class
class Rectangle: public Shape {
   public:
      int getArea() { 
         return (width * height); 
      };
      struct drift;
      void test(drift arg);
};

struct Rectangle::drift{
    size_t lhs_sv;
    double coeff;
    size_t rhs_sv;
    double order;
    };    
*/

#endif /* STOCHASTICMODEL_HPP_ */
