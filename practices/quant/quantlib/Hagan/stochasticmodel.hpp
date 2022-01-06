/*
 * stochasticmodel.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef STOCHASTICMODEL_HPP_
#define STOCHASTICMODEL_HPP_

#include <iostream>
#include <utility> // move
#include <memory> // unique_pntr

#include <ql/quantlib.hpp>
#include <ql/time/calendar.hpp>
#include <ql/pricingengine.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>
#include <Eigen/Dense> // Dense matrices

#include "model.hpp"
#include "optimizer.hpp"

// using namespace QuantLib;
namespace ql = QuantLib;

namespace StochasticModel{
   class HaganNF{
      public:
         HaganNF(boost::shared_ptr<ql::YieldTermStructure>& yieldCurve, const size_t& nFactor, Eigen::MatrixXd& corrMat);
         Eigen::MatrixXd getLowerMat() const {return lowerMat_;}
         Eigen::MatrixXd getCorrMat() const {return corrMat_;}
         void evolve(Eigen::MatrixXd& xn, const double& t, const Eigen::MatrixXd& x, const double& dt, const Eigen::MatrixXd& dw) const;
         Eigen::VectorXd evolve(const double& t, const Eigen::VectorXd& x, const double& dt, const Eigen::VectorXd& dw) const;
         std::shared_ptr<Eigen::MatrixXd>  computeInterestRate(const double& t_i, const double& dt,  const size_t& numPaths, const size_t& numSteps, std::vector<Eigen::MatrixXd>& x) const;
         double impliedVol(const ql::Period& today, const ql::Period& swaptionExpiry, const ql::Period& swaptionTenor, const double& tau=0.25, const ql::VolatilityType& type=ql::Normal) const;
         void calibrate(const std::shared_ptr<std::vector<ql::Period>>& swaptionExpiry, const std::shared_ptr<std::vector<ql::Period>>& swaptionTenor, const std::shared_ptr<Eigen::MatrixXd>& swaptionVolMat);
         void updateDZeta();
         // The followings are public variables due to calibration
         std::vector<Model::PolyFunc> alp_; // alpha (nFactor)
         std::vector<Model::ExpFunc> H_; // How to implement H_ in the framework? (nFactor)
         std::vector<std::vector<Model::PolyFunc>> dZeta_; // dZeta: alp^2. variable to compute integration of square of alphas. (nFactor, nFactor)

      private:
         size_t nFactor_;
         Eigen::MatrixXd corrMat_; // correlation matrices b.t.w. factors (nFactor, nFactor)
         Eigen::MatrixXd lowerMat_; // Lower part of the Cholesky decomposition of corrMat_ (nFactor, nFactor)
         boost::shared_ptr<ql::YieldTermStructure> yieldCurve_; // computes discount factor and forward rate
         struct HaganFunctor;

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
