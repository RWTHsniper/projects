#include "geometry.hpp"

namespace Geometry{
    void PolynomialCurve::computeFitting(){
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
    double PolynomialCurve::polyEval(double x) {
        double result = 0.0;
        for (int i = 0; i < coeffs_.size(); i++) {
            result += coeffs_[i] * pow(x, i);
        }
        return result;
    }
}
