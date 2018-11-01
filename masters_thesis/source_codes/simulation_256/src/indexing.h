/*
 * indexing.h
 *
 *  Created on: Oct 27, 2017
 *      Author: jaeyong
 */

#ifndef INDEXING_H_
#define INDEXING_H_

#include <iostream>     // std::cout
#include <stdio.h>
#include <stdlib.h>
#include <vector>

namespace helpers {

class indexing {
public:
	indexing(size_t *N,unsigned int dim);
	virtual ~indexing();
	// public methods
	size_t g(size_t I,size_t J);
	size_t g(size_t I,size_t J,size_t K);
	size_t g(size_t I,size_t J,size_t K,size_t L);
	size_t g(size_t I,size_t J,size_t K,size_t L,size_t M);
	size_t g(size_t I,size_t J,size_t K,size_t L,size_t M,size_t N);
	size_t g(size_t I,size_t J,size_t K,size_t L,size_t M,size_t N,size_t OO);


private:
	// Offering upto 6th order array.
	std::vector<size_t> n;
	//	unsigned int n[6];
};

} /* namespace helpers */

#endif /* INDEXING_H_ */
