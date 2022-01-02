#include "stochasticmodel.hpp"


void Rectangle::test(drift arg){
    std::cout << "hello" << std::endl;
    std::cout << "lhs_sv " << arg.lhs_sv << std::endl;
    std::cout << "coeff " << arg.coeff << std::endl;
    std::cout << "rhs_sv " << arg.rhs_sv << std::endl;
    std::cout << "order " << arg.order << std::endl;
}