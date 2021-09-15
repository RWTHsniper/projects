
#include "models.hpp"

GBM::GBM(double arg_r, double arg_sig) {
    num_sv = 1; // One price process
    r = arg_r;
    sig = arg_sig;
    std::cout << "Construct GBM " << r << " " << sig << std::endl;
};

GBM::~GBM() {
	// TODO Auto-generated destructor stub
};