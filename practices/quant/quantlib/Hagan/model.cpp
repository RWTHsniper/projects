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
    double PolyFunc::evaluate(double x) {
        double result = 0.0;
        for (int i = 0; i < coeffs_.size(); i++) {
            result += coeffs_[i] * pow(x, i);
        }
        return result;
    }
    // Evaluate the derivative of a polynomial.
    double PolyFunc::evalDeriv(double x) {
        double result = 0.0;
        for (int i = 1; i < coeffs_.size(); i++) {
            result += coeffs_[i] * i * pow(x, i-1);
        }
        return result;
    }

    void ExpFunc::computeFitting(){
        ExpFunc::check();
        for (size_t i=0; i<xVals_.size(); i++){
            buffer_[i] = evaluate(xVals_[i]);
        }

    };
    double ExpFunc::evaluate(double x) {
        double result = params_[0] * exp(params_[1] * x) + params_[2];
        return result;
    }
    double ExpFunc::evalDeriv(double x) {
        double result = params_[0] * params_[1] * exp(params_[1] * x);
        return result;
    }
}
