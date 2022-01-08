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
         HaganNF(boost::shared_ptr<ql::YieldTermStructure>& yieldCurve, const size_t& nFactor, Eigen::MatrixXd& corrMat, const size_t& alpOrder);
         Eigen::MatrixXd getLowerMat() const {return lowerMat_;}
         Eigen::MatrixXd getCorrMat() const {return corrMat_;}
         void evolve(Eigen::MatrixXd& xn, const double& t, const Eigen::MatrixXd& x, const double& dt, const Eigen::MatrixXd& dw) const;
         Eigen::VectorXd evolve(const double& t, const Eigen::VectorXd& x, const double& dt, const Eigen::VectorXd& dw) const;
         double computeNumeraire(const double& t, const Eigen::VectorXd& x) const;
         double computeDiscount(const double& t, const double& T, const Eigen::VectorXd& x) const;
         double computeForward(const double& t, const double& T1, const double& T2, const Eigen::VectorXd& x) const;
         std::shared_ptr<Eigen::MatrixXd>  computeInterestRate(const double& t_i, const double& dt,  const size_t& numPaths, const size_t& numSteps, const std::vector<Eigen::MatrixXd>& x) const;
         double impliedVol(const ql::Period& swaptionExpiry, const ql::Period& swaptionTenor, const double& tau=0.25, const ql::VolatilityType& type=ql::Normal) const;
         double impliedVolAnal(const ql::Period& swaptionExpiry, const ql::Period& swaptionTenor, const double& tau=0.25, const ql::VolatilityType& type=ql::Normal) const;
         void calibrate(const std::shared_ptr<std::vector<ql::Period>>& swaptionExpiry, const std::shared_ptr<std::vector<ql::Period>>& swaptionTenor, const std::shared_ptr<Eigen::MatrixXd>& swaptionVolMat);
         void updateDZeta();
         void getInfo();
         // The followings are public variables due to calibration
         std::vector<Model::PolyFunc> alp_; // alpha (nFactor)
         std::vector<Model::ExpFunc> H_; // How to implement H_ in the framework? (nFactor)
         std::vector<std::vector<Model::PolyFunc>> dZeta_; // dZeta: alp^2. variable to compute integration of square of alphas. (nFactor, nFactor)

      private:
         Eigen::VectorXd computeH(const double& t) const;
         Eigen::MatrixXd computeZeta(const double& t) const;
         size_t nFactor_;
         size_t alpOrder_;
         Eigen::MatrixXd corrMat_; // correlation matrices b.t.w. factors (nFactor, nFactor)
         Eigen::MatrixXd lowerMat_; // Lower part of the Cholesky decomposition of corrMat_ (nFactor, nFactor)
         boost::shared_ptr<ql::YieldTermStructure> yieldCurve_; // computes discount factor and forward rate
         struct HaganFunctor;

   };


   /*
   * Bachelier's Normal model for pricing swaptions
   */
   // inline double bachelierATMSwaption(double t, double T0, double tau, double sig, ql::Swap::Type type, std::vector<double>& discountFactors, double notional=1.0){
   inline double bachelierATMSwaption(const double& t, const double& T0, const double& tau, const double& sig, 
                                       Eigen::VectorXd discountFactors, ql::Swap::Type type=ql::Swap::Payer, const double& notional=1.0){
   /*
   Arguments
   - T0: expiration/reset date
   - t: current time
   - sig: implied volatility
   - type: type of a swaption (Payer or Receiver)
   - discountFactors: vector of discount factors at coupon payments
   - notional: notional amount
   */
      double payer_val{0.0};
      // Since the goal is to compute ATM swaption, R_swap(t) = K
      double sqrt_t = std::sqrt(T0-t); // sqrt of time to expiry
      double pdf_0{0.3989422804014327}; // pdf of Gaussian distribution at 0
      // payer_val = std::accumulate(discountFactors.begin(), discountFactors.end(), 0.0);
      payer_val = discountFactors.sum();
      payer_val *= notional * tau * sig * sqrt_t * pdf_0;

      switch (type) {
         case ql::Swap::Payer: return payer_val;
         // case ql::Swap::Receiver: return -payer_val;
         default: std::cout << "d"; exit(-1); // there are no applicable constant_expressions 
                                    // therefore default is executed
      }
   }

   struct BachiAtmSwaptionFunctor : Optimizer::Functor<double> {
      // members
      double npv_;
      double t_;
      double T0_;
      double tau_;
      Eigen::VectorXd discountFactors_;
      double initialSig_;
      // Simple constructor
        BachiAtmSwaptionFunctor(const double& npv, const double& t, const double& T0, const double& tau, const Eigen::VectorXd& discountFactors, const double& initialSig): 
            npv_(npv), t_(t), T0_(T0), tau_(tau), discountFactors_(discountFactors), initialSig_(initialSig), Optimizer::Functor<double>(1, 1) {};

         // Implementation of the objective function
         int operator () (const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {
            // std::cout << "current guess and output" << std::endl << z.transpose() << std::endl;
            fvec[0] = bachelierATMSwaption(t_, T0_, tau_, z[0], discountFactors_) - npv_;
            // std::cout << fvec.transpose() << std::endl;

            return 0;
        }
   };

   inline double computeIvol(const double& npv, const double& t, const double& T0, const double& tau, const Eigen::VectorXd& discountFactors, const double& initialSig){
      BachiAtmSwaptionFunctor functor(npv, t, T0, tau, discountFactors, initialSig);
      Eigen::NumericalDiff<BachiAtmSwaptionFunctor> numDiff(functor);
      Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BachiAtmSwaptionFunctor>,double> lm(numDiff);
      lm.parameters.maxfev = 1000;
      lm.parameters.xtol = 1.0e-10;
      Eigen::VectorXd z(1); z[0] = initialSig;
      int ret = lm.minimize(z);
      std::cout << "iter count: " << lm.iter << std::endl;
      std::cout << "return status: " << ret << std::endl; // status 2 is good
      Optimizer::LMReturnStatus(ret);
      // switch(ret){
      //    case Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation  : std::cout << "Too many function evaluations\n";   break;
      //    case Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall  : std::cout << "Relative reduction is too small\n";   break;
      //    case Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall  : std::cout << "Relative error is too small\n";   break;
      //    case Eigen::LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall  : std::cout << "Relative error and reduction are too small\n";   break;
      // }
      Eigen::VectorXd tmp(1);
      functor(z, tmp);
      std::cout << "Norm of final function: " << tmp.norm() << std::endl;
      std::cout << "L1-norm of final function: " << tmp.lpNorm<1>() << std::endl;
      std::cout << "zSolver: " << z.transpose() << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      return z[0];
    }
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
