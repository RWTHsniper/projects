#include "stochasticmodel.hpp"

namespace StochasticModel{
        HaganNF::HaganNF(boost::shared_ptr<ql::YieldTermStructure>& yieldCurve, const size_t& nFactor, Eigen::MatrixXd& corrMat, const size_t& alpOrder): yieldCurve_(yieldCurve),
            nFactor_(nFactor), corrMat_(corrMat), alpOrder_(alpOrder){
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
            // initialize alphas
            // Assume linear functions for alpha (volatility)
            Eigen::VectorXd polyCoeffs(alpOrder_+1); 
            if (alpOrder_==1) polyCoeffs << 0.0,1.0; // f(x) = x
            else if (alpOrder_ == 2) polyCoeffs << 0.0, 0.0, 1.0; // f(x) = x^2
            alp_.reserve(nFactor_);
            for (size_t i=0; i<nFactor_; i++){alp_.emplace_back(alpOrder_, polyCoeffs);}
            dZeta_.reserve(nFactor_);
            for (size_t i=0; i< nFactor_; i++){
                std::vector<Model::PolyFunc> tmp;
                tmp.reserve(nFactor_);
                for (size_t j=0; j< nFactor_; j++){
                    Model::PolyFunc elem = alp_[i] * alp_[j] * corrMat_(i, j);
                    tmp.emplace_back(elem); 
                }
                dZeta_.emplace_back(tmp);
            }
            H_.reserve(nFactor);
            double kappa = 1e-3;
            Eigen::VectorXd expParams(3); expParams << -1.0/kappa, -kappa, 1.0/kappa; // g(x) = (1-exp(-k*x))/k
            for (size_t i=0; i< nFactor_; i++){
                H_.emplace_back(expParams);
            }
            // Information to check
            // for (size_t i=0; i<nFactor_; i++){
            //    for (size_t j=0; j<nFactor_; j++)
            //       dZeta_[i][j].getInfo();
            // }
            // H_[0].getInfo();
            // H_[1].getInfo();
         }

    void HaganNF::evolve(Eigen::MatrixXd& xn, const double& t, const Eigen::MatrixXd& x, const double& dt, const Eigen::MatrixXd& dw) const{
        // Eigen::VectorXd xn(x); // state variables at next step (output)
        // Compute drift
        // LGM model does not have any drift
        // Compute volatility
        // Eigen::MatrixXd vol(xn.size());
        // std::cout << "Start stochastic" << std::endl;
        // std::cout << "rows and cols " << x.rows() << " " << x.cols() << std::endl;
        double vol = 0.0;
        double sqrtDt = std::sqrt(dt);
        for (size_t j=0; j< x.cols(); j++){
            for (size_t i=0; i< x.rows(); i++){
                // vol(i,j) = alp_[i].evaluate(t);
                vol = alp_[i].evaluate(t);
                xn(i,j) = x(i,j) + vol * dw(i,j) * sqrtDt;
            }
        }
        // for (size_t i=0; i< nFactor_; i++){
        //     vol[i] = alp_[i].evaluate(t);
        // }
        // for (size_t i=0; i< nFactor_; i++){
        //     xn[i] += vol[i] * dw[i] * std::sqrt(dt);
        // }
    }
    void HaganNF::updateDZeta(){
        /* Update dZeta_ after alp_ is updated
        *
        */
        for (size_t i=0; i< nFactor_; i++){
            for (size_t j=0; j< nFactor_; j++){
                dZeta_[i][j] = alp_[i] * alp_[j] * corrMat_(i, j);
            }
        }
    }

    Eigen::VectorXd HaganNF::evolve(const double& t, const Eigen::VectorXd& x, const double& dt, const Eigen::VectorXd& dw) const{
        Eigen::VectorXd xn(x); // state variables at next step (output)
        double sqrtDt = std::sqrt(dt);
        // Compute drift
        // LGM model does not have any drift
        // Compute volatility
        Eigen::VectorXd vol(nFactor_); 
        for (size_t i=0; i< nFactor_; i++){
            vol[i] = alp_[i].evaluate(t);
        }
        for (size_t i=0; i< nFactor_; i++){
            xn[i] += vol[i] * dw[i] * sqrtDt;
        }
        return xn;
    }
    std::shared_ptr<Eigen::MatrixXd> HaganNF::computeInterestRate(const double& t_i, const double& dt, const size_t& numPaths, const size_t& numSteps, std::vector<Eigen::MatrixXd>& x) const {
        double t = t_i + dt; // skip 0th step
        double T = t_i + numSteps * dt;
        std::shared_ptr<Eigen::MatrixXd> r = std::make_shared<Eigen::MatrixXd>(numPaths, numSteps+1);
        r->setZero();
        // Eigen::MatrixXd r(numPaths, numSteps+1); r.setZero();
        // temporary variables for H, H', and zeta
        Eigen::VectorXd H_t(nFactor_);
        Eigen::VectorXd H_t_deriv(nFactor_);
        Eigen::MatrixXd zeta_t(nFactor_, nFactor_);
        for (size_t step=1; step <= numSteps; step++){ // index for steps
            double fwdRate = yieldCurve_->forwardRate(0.0,6.0, ql::Continuous);
            for (size_t sv=0; sv<nFactor_; sv++){ // index for state variable
                H_t[sv] = H_[sv].evaluate(t);
                H_t_deriv[sv] = H_[sv].evalDeriv(t);
                for (size_t sv2=0; sv2<nFactor_; sv2++){ // index for state variable
                    zeta_t(sv, sv2) = dZeta_[sv][sv2].evalInt(t_i, t);
                }
            }
            for (size_t path=0; path < numPaths; path++){ // index for paths
                (*r)(path, step) = fwdRate;
                (*r)(path, step) += H_t_deriv.dot(x[step].col(path) + zeta_t * H_t);
            }
            // // dw = lowerMat_ * dWIndep[i];
            // dw = lowerMat_; dw *= dWIndep[i];
            // haganModel.evolve(x[i+1], t, x[i], dt, dw);
            // if (i==numSteps-1) saveData(source_dir+"output/dw.csv", dw);
            t+= dt;
        }
    return r;
    }
    double HaganNF::impliedVol(const ql::Period& today, const ql::Period& swaptionExpiry, const ql::Period& swaptionTenor, const double& tau, const ql::VolatilityType& type) const { // type: ql::Normal
        double res = 0.0;
        double t = static_cast<double>(ql::years(today));
        Eigen::MatrixXd zeta_t(nFactor_, nFactor_); zeta_t.setZero();
        size_t N = static_cast<size_t>(ql::years(swaptionTenor) / tau); // number of payments for a swap
        double t0 = ql::years(swaptionExpiry); // exercise date
        double ti = t0;
        double tn = t0 + ql::years(swaptionTenor);
        Eigen::VectorXd Pi(N+1); // t0,t1,...,tN
        Eigen::VectorXd dHi(nFactor_);
        Eigen::VectorXd Htot(nFactor_); Htot.setZero();
        // compute discount factors
        for (size_t i=0; i<N+1; i++){ // t0,t1,...,tN
            Pi[i] = yieldCurve_->discount(ti); // equivalent to Di
            ti += tau;
        }
        // compute annuity
        double annuity = 0.0;
        for (size_t i=1; i<=N; i++){ // t1,...,tN
            annuity += tau * Pi[i];
        }
        double St = (Pi[0] - Pi[N]) / annuity;
        if (type == ql::Normal){
            // std::cout << t0 << " " << ti << " " << tn << " " << std::endl << Pi << std::endl;
            // std::cout << N << std::endl;
            // compute Zeta
            for (size_t i=0; i<nFactor_; i++)
                for (size_t j=0; j<nFactor_; j++){
                    zeta_t(i,j) = dZeta_[i][j].evalInt(0.0, t0);
                }
            // compute Htot
            for (size_t j=0; j<nFactor_; j++){
                ti = t0;
                Htot[j] = Pi[N] * (H_[j].evaluate(tn) - H_[j].evaluate(t));
                for (size_t i=1; i<=N; i++){
                    ti = t0 + tau;
                    Htot[j] += St * tau * Pi[i] * (H_[j].evaluate(ti) - H_[j].evaluate(t));
                }
            }
            Htot /= annuity;
            // sigN^2 * tex = Htot.dot(Zeta*Htot)
            res = std::sqrt(Htot.dot(zeta_t*Htot) / t0);
        }
        else{
            std::cout << "Volatility type " << type << " is not supported" << std::endl;
        }
        return res;
    }
    struct HaganNF::HaganFunctor : Optimizer::Functor<double> {
        // members
        HaganNF* this_;
        std::shared_ptr<std::vector<ql::Period>> swaptionExpiry_;
        std::shared_ptr<std::vector<ql::Period>> swaptionTenor_;
        std::shared_ptr<Eigen::MatrixXd> swaptionVolMat_; 
        // Simple constructor
        HaganFunctor(HaganNF* thisObj, 
                    const std::shared_ptr<std::vector<ql::Period>>& swaptionExpiry, 
                    const std::shared_ptr<std::vector<ql::Period>>& swaptionTenor,
                    const std::shared_ptr<Eigen::MatrixXd>& swaptionVolMat): 
            this_(thisObj), swaptionExpiry_(swaptionExpiry), swaptionTenor_(swaptionTenor), swaptionVolMat_(swaptionVolMat), 
            Optimizer::Functor<double>(swaptionExpiry->size()*swaptionTenor->size(),swaptionVolMat->rows()*swaptionVolMat->cols()) {}

        // Implementation of the objective function
        int operator () (const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {
            ql::Period today(0, ql::Years); // Current time for evaluation
            ql::VolatilityType type = ql::Normal;
            double tau = 0.25;
            /*
            * Evaluate the fitness to ATM swaption matrix.
            * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
            * of squared terms. The algorithm takes this into account: do not do it yourself.
            * In other words: objFun = sum(fvec(i)^2)
            */
           // Setup the parameters in LGM
            double penalty = 0.0; // penalty function
            size_t count = 0;
            for (size_t i=0; i<this_->nFactor_; i++){
                size_t order = this_->alp_[i].getOrder();
                // auto coeffs = this_->alp_[i].getCoeffs();
                Eigen::VectorXd buff(order+1);
                for (size_t j=0; j<= order; j++){ // coefficients from 0th to the last order
                    buff[j] = z[count];
                    count++;
                }
                this_->alp_[i].setCoeffs(buff);
            }
            for (size_t i=0; i<this_->nFactor_; i++){
                Eigen::VectorXd buff(3);
                auto expParams = this_->H_[i].getParams(); // a*exp(b) + c: [a, b, c]: [-1/ka, -ka, 1/ka]
                double kappa = z[count];
                buff << -1/kappa, -kappa, 1/kappa;
                this_->H_[i].setParams(buff);
                count++;
                if (kappa <= 0.0){
                    // penalty += std::pow(kappa, 2);
                }
                else if (kappa > 10){
                    penalty += std::pow(kappa-10.0, 2);
                }
            }
            this_->updateDZeta();
            // std::cout << "current sol" << std::endl << z.transpose() << std::endl;
            for (size_t i=0; i< swaptionExpiry_->size(); i++){
                for (size_t j=0; j< swaptionTenor_->size(); j++){
                    size_t ind = i*swaptionTenor_->size() + j;
                    // Use z (xvalues) to update model parameters
                    double iVol = this_->impliedVol(today, (*swaptionExpiry_)[i], (*swaptionTenor_)[j], tau, type);
                    fvec(ind) = (*swaptionVolMat_)(i,j) - iVol;
                    // std::cout << (*swaptionVolMat_)(i,j) - iVol << " ";
                    fvec(ind) += penalty / fvec.size();
                }
                // std::cout << std::endl;                
            }
            return 0;
        }
    };

    void HaganNF::calibrate(const std::shared_ptr<std::vector<ql::Period>>& swaptionExpiry, 
                            const std::shared_ptr<std::vector<ql::Period>>& swaptionTenor,
                            const std::shared_ptr<Eigen::MatrixXd>& swaptionVolMat){
        /* The parameters to calibrate are
        * alp_[0], alp_[1], ..., alp_[nFactor_-1]: coeff[0], coeff[1]
        * H_[0], ..., H_[nFactor-1]: k 
        * In total, nFactor * 2 + nFactor * 1 parameters should be updated
        */
        size_t numParams = 0;
        for (size_t i=0; i< alp_.size(); i++) numParams += alp_[i].getOrder() + 1; // number of polynomial coefficients
        numParams += nFactor_; // one param for H_
        // size_t numParams = nFactor_ * (2); // alp_ and H_
        // Initialize the vector of initial values
        Eigen::VectorXd zInit(numParams);
        size_t count = 0;
        for (size_t i=0; i<nFactor_; i++){
            auto coeffs = alp_[i].getCoeffs();
            size_t order = alp_[i].getOrder();
            for (size_t j=0; j<= order; j++){ // coefficients from 0th to the last order
                zInit[count] = coeffs[j];
                count++;
            }
        }
        for (size_t i=0; i<nFactor_; i++){
            auto expParams = H_[i].getParams(); // a*exp(b) + c: [a, b, c]
            double kappa = -expParams[1]; 
            zInit[count] = kappa;
            count++;
        }
        std::cout << "Initial guess for LGM " << std::endl << zInit.transpose() << std::endl;
 
        HaganFunctor functor(this, swaptionExpiry, swaptionTenor, swaptionVolMat);
        Eigen::NumericalDiff<HaganFunctor> numDiff(functor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<HaganFunctor>,double> lm(numDiff);
        lm.parameters.maxfev = 1000;
        lm.parameters.xtol = 1.0e-10;
        std::cout << "max fun eval: " << lm.parameters.maxfev << std::endl;
        std::cout << "x tol: " << lm.parameters.xtol << std::endl;

        Eigen::VectorXd z = zInit;
        int ret = lm.minimize(z);
        std::cout << "iter count: " << lm.iter << std::endl;
        std::cout << "return status: " << ret << std::endl; // status 2 is good
        switch(ret){
            case Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation  : std::cout << "Too many function evaluations\n";   break;
            case Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall  : std::cout << "Relative reduction is too small\n";   break;
            case Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall  : std::cout << "Relative error is too small\n";   break;
            case Eigen::LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall  : std::cout << "Relative error and reduction are too small\n";   break;
        }
        Eigen::VectorXd tmp(swaptionVolMat->rows()*swaptionVolMat->cols());
        functor(z, tmp);
        std::cout << "Norm of final function: " << tmp.norm() << std::endl;
        std::cout << "L1-norm of final function: " << tmp.lpNorm<1>() << std::endl;
        std::cout << "zSolver: " << z.transpose() << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }

}

/*
namespace LevenbergMarquardtSpace {
    enum Status {
        NotStarted = -2,
        Running = -1,
        ImproperInputParameters = 0,
        RelativeReductionTooSmall = 1,
        RelativeErrorTooSmall = 2,
        RelativeErrorAndReductionTooSmall = 3,
        CosinusTooSmall = 4,
        TooManyFunctionEvaluation = 5,
        FtolTooSmall = 6,
        XtolTooSmall = 7,
        GtolTooSmall = 8,
        UserAsked = 9
    };
}
*/