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
                                    Eigen::VectorXd polyCoeffs(2); polyCoeffs << 0.0,1.0; // f(x) = x
                                    // initialize alphas
                                    alp.reserve(nFactor_);
                                    size_t order = 1;
                                    for (size_t i=0; i<nFactor_; i++){alp.emplace_back(order, polyCoeffs);}
                                    zeta_.reserve(nFactor_);
                                    for (size_t i=0; i< nFactor_; i++){
                                       std::vector<Model::PolyFunc> tmp;
                                       tmp.reserve(nFactor_);
                                       for (size_t j=0; j< nFactor_; j++){
                                          Model::PolyFunc elem = alp[i] * alp[j] * corrMat_(i, j);
                                          tmp.emplace_back(elem); 
                                       }
                                       zeta_.emplace_back(tmp);
                                    }
                                    H_.reserve(nFactor);
                                    double kappa = 1.0;
                                    Eigen::VectorXd expParams(3); expParams << -1.0/kappa, -kappa, 1.0/kappa; // g(x) = (1-exp(-k*x))/k
                                    for (size_t i=0; i< nFactor_; i++){
                                       H_.emplace_back(expParams);
                                    }
                                    // Information to check
                                    // for (size_t i=0; i<nFactor_; i++){
                                    //    for (size_t j=0; j<nFactor_; j++)
                                    //       zeta_[i][j].getInfo();
                                    // }
                                    // H_[0].getInfo();
                                    // H_[1].getInfo();
         }
         Eigen::MatrixXd getLowerMat() const {return lowerMat_;}
         Eigen::MatrixXd getCorrMat() const {return corrMat_;}
         Eigen::VectorXd evolve(const double& t, const Eigen::VectorXd& x, const double& dt, const Eigen::VectorXd& dw) const;

      private:
         size_t nFactor_;
         std::vector<Model::PolyFunc> alp; // alpha (nFactor)
         std::vector<Model::ExpFunc> H_; // How to implement H_ in the framework? (nFactor)
         std::vector<std::vector<Model::PolyFunc>> zeta_; // integration fo square of alphas. (nFactor, nFactor)
         Eigen::MatrixXd corrMat_; // correlation matrices b.t.w. factors (nFactor, nFactor)
         Eigen::MatrixXd lowerMat_; // Lower part of the Cholesky decomposition of corrMat_ (nFactor, nFactor)
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
