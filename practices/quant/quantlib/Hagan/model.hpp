/*
 * models.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <iostream>
#include <Eigen/Dense> // Dense matrices
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


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
            double evalDeriv(const double& x, const size_t& order) const;
            double evalInt(const double& x_i, const double& x_f) const;
            void getInfo() const;
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
            void setxVals(const Eigen::VectorXd& xVals){xVals_ = xVals;}
            void setyVals(const Eigen::VectorXd& yVals){yVals_ = yVals;}
            void setParams(const Eigen::VectorXd& params){assert(params_.size() == 3); params_ = params;}
            double evaluate(const double& x);
            double evalDeriv(const double& x, const size_t& order);
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
    
// // https://stackoverflow.com/questions/18509228/how-to-use-the-eigen-unsupported-levenberg-marquardt-implementation
// // Generic functor
// // See http://eigen.tuxfamily.org/index.php?title=Functors
// // C++ version of a function pointer that stores meta-data about the function
// template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
// struct Functor
// {

//   // Information that tells the caller the numeric type (eg. double) and size (input / output dim)
//   typedef _Scalar Scalar;
//   enum { // Required by numerical differentiation module
//       InputsAtCompileTime = NX,
//       ValuesAtCompileTime = NY
//   };

//   // Tell the caller the matrix sizes associated with the input, output, and jacobian
//   typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
//   typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
//   typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

//   // Local copy of the number of inputs
//   int m_inputs, m_values;

//   // Two constructors:
//   Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
//   Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

//   // Get methods for users to determine function input and output dimensions
//   int inputs() const { return m_inputs; }
//   int values() const { return m_values; }

// };



// // https://en.wikipedia.org/wiki/Test_functions_for_optimization
// // Booth Function
// // Implement f(x,y) = (x + 2*y -7)^2 + (2*x + y - 5)^2
// struct BoothFunctor : Functor<double>
// {
//   // Simple constructor
//   BoothFunctor(): Functor<double>(2,2) {}

//   // Implementation of the objective function
//   int operator()(const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {
//     double x = z(0);   double y = z(1);
//     /*
//      * Evaluate the Booth function.
//      * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
//      * of squared terms. The algorithm takes this into account: do not do it yourself.
//      * In other words: objFun = sum(fvec(i)^2)
//      */
//     fvec(0) = x + 2*y - 7;
//     fvec(1) = 2*x + y - 5;
//     return 0;
//   }
// };


// // https://en.wikipedia.org/wiki/Test_functions_for_optimization
// // Himmelblau's Function
// // Implement f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
// struct HimmelblauFunctor : Functor<double>
// {
//   // Simple constructor
//   HimmelblauFunctor(): Functor<double>(2,2) {}

//   // Implementation of the objective function
//   int operator()(const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {
//     double x = z(0);   double y = z(1);
//     /*
//      * Evaluate Himmelblau's function.
//      * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
//      * of squared terms. The algorithm takes this into account: do not do it yourself.
//      * In other words: objFun = sum(fvec(i)^2)
//      */
//     fvec(0) = x * x + y - 11;
//     fvec(1) = x + y * y - 7;
//     return 0;
//   }
// };


// void testBoothFun() {
//   std::cout << "Testing the Booth function..." << std::endl;
//   Eigen::VectorXd zInit(2); zInit << 1.87, 2.032;
//   std::cout << "zInit: " << zInit.transpose() << std::endl;
//   Eigen::VectorXd zSoln(2); zSoln << 1.0, 3.0;
//   std::cout << "zSoln: " << zSoln.transpose() << std::endl;

//   BoothFunctor functor;
//   Eigen::NumericalDiff<BoothFunctor> numDiff(functor);
//   Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BoothFunctor>,double> lm(numDiff);
//   lm.parameters.maxfev = 1000;
//   lm.parameters.xtol = 1.0e-10;
//   std::cout << "max fun eval: " << lm.parameters.maxfev << std::endl;
//   std::cout << "x tol: " << lm.parameters.xtol << std::endl;

//   Eigen::VectorXd z = zInit;
//   int ret = lm.minimize(z);
//   std::cout << "iter count: " << lm.iter << std::endl;
//   std::cout << "return status: " << ret << std::endl;
//   std::cout << "zSolver: " << z.transpose() << std::endl;
//   std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
// }


}


#endif /* MODEL_HPP_ */
