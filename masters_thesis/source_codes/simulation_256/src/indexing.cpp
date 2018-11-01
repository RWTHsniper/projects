/*
 * indexing.cpp
 *
 *  Created on: Oct 27, 2017
 *      Author: jaeyong
 */

#include "indexing.h"

namespace helpers {

indexing::indexing(size_t* N,unsigned int dim) {
	// TODO Auto-generated constructor stub

	n.resize(dim,0);
	if (dim > 7){
		std::cout <<"indexing dimension exceeding 7 not supported. Input dimension is "<<dim<<std::endl;
		abort();
	}
	for (size_t i=0;i<dim;i++) n[i]=N[i];
}

size_t indexing::g(size_t I,size_t J){
	if (n.size() != 2)	std::cout << "Wrong indexing!! n.size() and num indx "<<n.size()<<" "<<2<< std::endl;
	size_t outp = I*n[1]+J;
	return outp;
}
size_t indexing::g(size_t I,size_t J,size_t K){
	if (n.size() != 3)	std::cout << "Wrong indexing!! n.size() and num indx "<<n.size()<<" "<<3<< std::endl;
	size_t outp = I*n[1]*n[2]+J*n[2]+K;
	return outp;
}
size_t indexing::g(size_t I,size_t J,size_t K,size_t L){
	if (n.size() != 4)	std::cout << "Wrong indexing!! n.size() and num indx "<<n.size()<<" "<<4<< std::endl;
	size_t outp = I*n[1]*n[2]*n[3]+J*n[2]*n[3]+K*n[3]+L;
	return outp;
}
size_t indexing::g(size_t I,size_t J,size_t K,size_t L,size_t M){
	if (n.size() != 5)	std::cout << "Wrong indexing!! n.size() and num indx "<<n.size()<<" "<<5<< std::endl;
	size_t outp = I*n[1]*n[2]*n[3]*n[4]+J*n[2]*n[3]*n[4]+K*n[3]*n[4]+L*n[4]+M;
	return outp;
}
size_t indexing::g(size_t I,size_t J,size_t K,size_t L,size_t M,size_t N){
	if (n.size() != 6)	std::cout << "Wrong indexing!! n.size() and num indx "<<n.size()<<" "<<6<< std::endl;
	size_t outp = I*n[1]*n[2]*n[3]*n[4]*n[5]+J*n[2]*n[3]*n[4]*n[5]+K*n[3]*n[4]*n[5]+L*n[4]*n[5]+M*n[5]+N;
	return outp;
}
size_t indexing::g(size_t I,size_t J,size_t K,size_t L,size_t M,size_t N,size_t OO){
	if (n.size() != 7)	std::cout << "Wrong indexing!! n.size() and num indx "<<n.size()<<" "<<7<< std::endl;
	size_t outp = I*n[1]*n[2]*n[3]*n[4]*n[5]*n[6]+J*n[2]*n[3]*n[4]*n[5]*n[6]+K*n[3]*n[4]*n[5]*n[6]+L*n[4]*n[5]*n[6]+M*n[5]*n[6]+N*n[6]+OO;
	return outp;
}

indexing::~indexing() {
	// TODO Auto-generated destructor stub

}

} /* namespace helpers */
