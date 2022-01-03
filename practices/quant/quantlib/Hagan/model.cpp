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

    PolyFunc PolyFunc::getIndefIntegral() const { // get indefinite integral
        
    } 


    void ExpFunc::computeFitting(){
        ExpFunc::check();
        for (size_t i=0; i<xVals_.size(); i++){
            buffer_[i] = evaluate(xVals_[i]);
        }

    };
    double ExpFunc::evaluate(const double& x) {
        double result = params_[0] * exp(params_[1] * x) + params_[2];
        return result;
    }
    double ExpFunc::evalDeriv(const double& x, const size_t& order) {
        assert(order > 0);
        double result = pow(params_[1], order) * params_[0] * exp(params_[1] * x);
        return result;
    }
}
