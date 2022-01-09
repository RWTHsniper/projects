#include "model.hpp"

namespace Model{
    void PolyFunc::computeFitting(){
        PolyFunc::check();
        for (size_t i = 0; i < xVals_.size(); i++) {
            A_(i, 0) = 1.0;
        }
        for (size_t j = 0; j < xVals_.size(); j++) {
            for (int i = 0; i < order_; i++) {
            A_(j, i + 1) = A_(j, i) * xVals_(j);
            }
        }
        auto& Q = A_.householderQr();
        coeffs_ = Q.solve(yVals_);
    };

    // Evaluate a polynomial.
    double PolyFunc::evaluate(const double& x) const {
        double result = 0.0;
        for (int i = 0; i < coeffs_.size(); i++) {
            result += coeffs_[i] * pow(x, i);
        }
        return result;
    }
    // Evaluate the derivative of a polynomial.
    double PolyFunc::evalDeriv(const double& x, const size_t& order) const {
        assert(order > 0);
        double result = 0.0;
        for (size_t i = order; i < coeffs_.size(); i++) { // position in coefficient vector
            double tmp = 1.0; // multiplication of powers
            for (size_t j = 0; j < order; j++){
                tmp *= i-j;
            }
            result += coeffs_[i] * tmp * pow(x, i-order);
        }
        return result;
    }
    double PolyFunc::evalInt(const double& x_i, const double& x_f) const {
        double result = 0.0;
        for (size_t n = 0; n < coeffs_.size(); n++) { // position and order for each coefficient
            result += coeffs_[n] / (n+1) * (pow(x_f, n+1) - pow(x_i, n+1));
        }
        return result;
    }

    void PolyFunc::getInfo() const {
        std::cout << "Information of a PolyFunc " << this << std::endl;
        std::cout << "order " << std::endl << order_ << std::endl;
        std::cout << "coeffs " << std::endl << coeffs_.transpose() << std::endl;
        std::cout << "expression" << std::endl;
        for (int i= order_; i>=0; i--){
            if (i != order_) {std::cout  << " + ";}
            if (i>0){
                std::cout << coeffs_[i] << "*x^(" << i << ")";
            }
            else{
                std::cout << coeffs_[i] << std::endl;
            }
        }
    }

    // PolyFunc PolyFunc::getIndefIntegral() const { // get indefinite integral
        
    // } 


    void ExpFunc::computeFitting(){
        ExpFunc::check();
        for (size_t i=0; i<xVals_.size(); i++){
            buffer_[i] = evaluate(xVals_[i]);
        }

    };
    double ExpFunc::evaluate(const double& x) const {
        double result = params_[0] * exp(params_[1] * x) + params_[2];
        return result;
    }
    double ExpFunc::evalDeriv(const double& x, const size_t& order) const {
        assert(order > 0);
        double result = pow(params_[1], order) * params_[0] * exp(params_[1] * x);
        return result;
    }
    double ExpFunc::evalInt(const double& x_i, const double& x_f) const {
        double result = params_[0] / params_[1] * (std::exp(params_[1] * x_f) - std::exp(params_[1] * x_f)) + params_[2]*(x_f - x_i);
        return result;
    }
    void ExpFunc::getInfo() const {
        std::cout << "Information of a ExpFunc " << this << std::endl;
        std::cout << params_[0] << "*Exp[" << params_[1] <<"*x] + " << params_[2] << std::endl;
    }

    void testModule(){
        double tol{1e-8};
        size_t order = 0;
        Eigen::VectorXd coeffs(order+1); coeffs << 2.5;
        PolyFunc poly0(order, coeffs);
        assert(tol >= abs(poly0.evaluate(5) - 2.5)); // should be less than a tolerance
        order = 1; 
        coeffs.resize(order+1); coeffs << 3, 5; // 5*x + 3
        PolyFunc poly1(order, coeffs);
        assert(tol >= abs(poly1.evaluate(10) - 53)); // should be less than a tolerance
        assert(tol >= abs(poly1.evalInt(0, 10) - 280)); // should be less than a tolerance
        assert(tol >= abs(poly1.evalDeriv(-20) - 5)); // should be less than a tolerance

        Eigen::VectorXd xVals(5);
        Eigen::VectorXd yVals(5);
        xVals << 1,2,3,4,5;
        yVals << 1,4,9,16,25;
        Model::PolyFunc polyCurve(xVals, yVals, 2);
        for (size_t i=0; i<xVals.size(); i++){
            std::cout << polyCurve.evaluate(xVals[i]) << " " << yVals[i] << std::endl;
            assert(tol >= std::abs(polyCurve.evaluate(xVals[i]) - yVals[i])); // should be less than a tolerance
        }
        std::cout << "+1" << std::endl;
        polyCurve += 1;
        for (size_t i=0; i<xVals.size(); i++){
            std::cout << polyCurve.evaluate(xVals[i]) << " " << yVals[i]+1 << std::endl;
            assert(tol >= std::abs(polyCurve.evaluate(xVals[i]) - (yVals[i]+1))); // should be less than a tolerance
        }
        std::cout << "Integration of x^2 + 1" << std::endl; // 1/3*x^3 + x 
        for (size_t i=0; i<xVals.size(); i++){
            std::cout << (1/3.0 * pow(xVals[i],3) + xVals[i]) << " " << polyCurve.evalInt(0.0, xVals[i]) << std::endl;
            assert(tol >= std::abs(polyCurve.evalInt(0.0, xVals[i]) - (1/3.0 * pow(xVals[i],3) + xVals[i])));
        }
        Model::PolyFunc polyCurve2(polyCurve);
        std::cout << "add1 " << &(polyCurve) << std::endl;
        std::cout << "add2 " << &(polyCurve2) << std::endl;
        Model::PolyFunc polyCurve3 = polyCurve + polyCurve2; // (x^2+1)*2
        std::cout << "add3 " << &(polyCurve3) << std::endl;
        for (size_t i=0; i<xVals.size(); i++){
            std::cout << polyCurve3.evaluate(xVals[i]) << " " << 2.0*(pow(xVals[i],2)+1) << std::endl;
            assert(tol >= std::abs(polyCurve3.evaluate(xVals[i]) - 2.0*(pow(xVals[i],2)+1)));
        }
        Model::PolyFunc polyCurve4 = polyCurve*polyCurve2; // (x^2+1)^2
        for (size_t i=0; i<xVals.size(); i++){
            std::cout << polyCurve4.evaluate(xVals[i]) << " " << pow(pow(xVals[i],2)+1,2) << std::endl;
            assert(tol >= std::abs(polyCurve4.evaluate(xVals[i]) - pow(pow(xVals[i],2)+1,2)));
        }

        Eigen::VectorXd params(3);
        params << 1.0, 1.0, 1.0;
        Eigen::VectorXd yExpVals(5); // exp(x) + 1
        yExpVals << exp(1)+1,exp(2)+1,exp(3)+1,exp(4)+1,exp(5)+1;
        Model::ExpFunc expCurve(xVals, yExpVals, params);
        for (size_t i=0; i<xVals.size(); i++){
            std::cout << expCurve.evaluate(xVals[i]) << " " << yExpVals[i] << std::endl;
            assert(tol >= std::abs(expCurve.evaluate(xVals[i]) - yExpVals[i])); // should be less than a tolerance
        }


        std::cout << "Test MODEL.HPP is complete" << std::endl;
    }

}
