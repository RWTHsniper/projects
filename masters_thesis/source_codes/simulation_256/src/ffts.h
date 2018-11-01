/*
 * ffts.h
 *
 *  Created on: Oct 30, 2017
 *      Author: jaeyong
 */

#ifndef FFTS_H_
#define FFTS_H_

#include <iostream>     // std::cout
#include <complex>
#include "fftw3.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

namespace fft{

void fftw(size_t* N,size_t dim,complex<double>* inp, complex<double>* outp,unsigned int dir,unsigned int nprocs);
void dft_coeffs(complex<double>* coefs,unsigned int order, size_t* N,double* L);	// build up coefficient array for FFT based differentiation.

}



#endif /* FFTS_H_ */
