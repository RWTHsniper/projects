#ifndef MODELS_HPP_
#define MODELS_HPP_

#include <iostream>
#include <vector>


class Models {
protected:
    size_t num_sv; // number of state variables for each model
public:
	virtual ~Models(){}; // empty destructor for a virtual class
};


class GBM: public Models {
// class GBM {
private:
	double r; // risk-free rate
	double sig; // volatility

public:
	GBM(double arg_r,  double arg_sig);
	virtual ~GBM();
};


#endif /* MODELS_HPP_ */