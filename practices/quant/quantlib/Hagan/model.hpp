/*
 * models.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <iostream>
#include <cassert> // assert
#include <Eigen/Dense> // Dense matrices


namespace Model{
    // Class definitions
    class PolyFunc{ // polynomial function
        /*
        https://gist.github.com/ksjgh/4d5050d0e9afc5fdb5908734335138d0
        */
        public:
            PolyFunc(const Eigen::VectorXd& xVals, const Eigen::VectorXd& yVals, const size_t& order):
                                                    xVals_(xVals), yVals_(yVals), order_(order){
                                                    A_.resize(xVals_.size(), order_ + 1);
                                                    coeffs_.resize(order_ + 1);
                                                    computeFitting();}
            PolyFunc(const size_t& order, const Eigen::VectorXd& coeffs): order_(order), coeffs_(coeffs){
                                                    xVals_.resize(0);
                                                    yVals_.resize(0);
            };
            void setCoeffs(const Eigen::VectorXd& coeffs){
                coeffs_ = coeffs;} // copy assignment            
            void setOrder(const size_t& order){order_ = order;}
            void setxVals(const Eigen::VectorXd& xVals){xVals_ = xVals;}
            void setyVals(const Eigen::VectorXd& yVals){yVals_ = yVals;}
            double evaluate(const double& x) const;
            double evalDeriv(const double& x, const size_t& order = 1) const;
            double evalInt(const double& x_i, const double& x_f) const;
            void getInfo() const;
            size_t getOrder() const {return order_;};
            Eigen::VectorXd getCoeffs(){return coeffs_;};
            PolyFunc getIndefIntegral() const; // get indefinite integral
            template <typename T>
            PolyFunc& operator+=(const T& other){
                this->coeffs_[0] += static_cast<double>(other);
                return *this;
            }
            template <typename T>
            PolyFunc& operator-=(const T& other){
                this->coeffs_[0] -= static_cast<double>(other);
                return *this;
            }
            PolyFunc operator+(const PolyFunc& other) const {
                PolyFunc res(*this); // copy constructor
                res.order_ = std::max(this->order_, other.order_);
                res.coeffs_.resize(res.order_+1);
                res.coeffs_.setZero();
                for(size_t i=0; i< this->coeffs_.size(); i++) res.coeffs_[i] += this->coeffs_[i];
                for(size_t i=0; i< other.coeffs_.size(); i++) res.coeffs_[i] += other.coeffs_[i];
                return res;
            }
            template <typename T>
            PolyFunc operator*(const T& other) const {
                PolyFunc res(*this); // copy constructor
                // for (size_t i=0; i<res.coeffs_.size(); i++){res.coeffs_[i] *= other;};
                res.coeffs_ *= other; // multiply by the factor
                return res;
            }
            // multiplication b.t.w. two polynomials
            PolyFunc operator*(const PolyFunc& other) const {
                PolyFunc res(*this); // copy constructor
                // To avoid any confusion, initialize xVals_ and yVals_ as empty vectors
                res.xVals_.resize(0);
                res.yVals_.resize(0);
                res.order_ = this->order_ + other.order_; // increas in the order
                res.coeffs_.resize(res.order_ + 1);
                res.coeffs_.setZero();
                // for (size_t i=0; i<this->coeffs_.size(); i++){
                //     for (size_t j=0; j<other.coeffs_.size(); j++){
                for (size_t i=0; i<this->order_+1; i++){
                    for (size_t j=0; j<other.order_+1; j++){
                        res.coeffs_[i+j] += this->coeffs_[i] * other.coeffs_[j];
                    }
                }
                return res;
            }
        private:
            size_t order_;
            Eigen::VectorXd xVals_;
            Eigen::VectorXd yVals_;
            Eigen::VectorXd coeffs_; // polynomial coefficients [0th, 1th, 2th, ...]
            Eigen::MatrixXd A_;
            void computeFitting();
            void check(){
                assert(xVals_.size() == yVals_.size());
                assert(order_ >= 1 && order_ <= xVals_.size() - 1);
            }
    };

    class ExpFunc{ // exponential function
        /*
        y(x) = a*exp(b*x) + c
        */
        public:
            ExpFunc(const Eigen::VectorXd& xVals, const Eigen::VectorXd& yVals, const Eigen::VectorXd& params):
                                                    xVals_(xVals), yVals_(yVals), params_(params){check(); buffer_.resize(xVals_.size());} // copy constructor
            ExpFunc(const Eigen::VectorXd& params): params_(params){
                                                    xVals_.resize(0);
                                                    yVals_.resize(0);
            };
            void setxVals(const Eigen::VectorXd& xVals){xVals_ = xVals;}
            void setyVals(const Eigen::VectorXd& yVals){yVals_ = yVals;}
            void setParams(const Eigen::VectorXd& params){assert(params_.size() == 3); params_ = params;}
            double evaluate(const double& x) const;
            double evalDeriv(const double& x, const size_t& order = 1) const;
            double evalInt(const double& x_i, const double& x_f) const;
            void getInfo() const;
            Eigen::VectorXd getParams() const {return params_;};

        private:
            Eigen::VectorXd xVals_;
            Eigen::VectorXd yVals_;
            Eigen::VectorXd params_; // [a, b, c]
            Eigen::VectorXd buffer_; // This buffer is size of yVals_.size() and contains intermediate values in fitting
            void computeFitting();
            void check(){
                assert(xVals_.size() == yVals_.size());
                assert(params_.size() == 3);
            }
    };
    
    void testModule();




}


#endif /* MODEL_HPP_ */
