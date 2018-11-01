/*
 * ffts.cpp
 *
 *  Created on: Nov 2, 2017
 *      Author: jaeyong
 */

#include "ffts.h"

namespace fft{

void fftw(size_t* N,size_t dim,complex<double>* inp, complex<double>* outp,unsigned int dir,unsigned int nprocs){

	if (nprocs >1){
		if (fftw_init_threads() == 0){cout<<"Error in initializing threads"<<endl;abort();}
		omp_set_dynamic(nprocs);
		omp_set_num_threads(nprocs);
		fftw_plan_with_nthreads(nprocs);
//		cout << "Running FFTW with "<<nprocs<<endl;
	}
	if (dir == 1){
		if (dim ==1){
			fftw_plan plan=fftw_plan_dft_1d(N[0],reinterpret_cast<fftw_complex*>(&inp[0]),reinterpret_cast<fftw_complex*>(&outp[0]),FFTW_FORWARD,FFTW_ESTIMATE);
			fftw_execute(plan);
		}
		else if (dim==2){
			fftw_plan plan=fftw_plan_dft_2d(N[0],N[1],reinterpret_cast<fftw_complex*>(&inp[0]),reinterpret_cast<fftw_complex*>(&outp[0]),FFTW_FORWARD,FFTW_ESTIMATE);
			fftw_execute(plan);
		}
		else if (dim==3){
			fftw_plan plan=fftw_plan_dft_3d(N[0],N[1],N[2],reinterpret_cast<fftw_complex*>(&inp[0]),reinterpret_cast<fftw_complex*>(&outp[0]),FFTW_FORWARD,FFTW_ESTIMATE);
			fftw_execute(plan);
		}
	}
	else if (dir == 0)
	{
		if (dim ==1){
			fftw_plan plan=fftw_plan_dft_1d(N[0],reinterpret_cast<fftw_complex*>(&inp[0]),reinterpret_cast<fftw_complex*>(&outp[0]),FFTW_BACKWARD,FFTW_ESTIMATE);
			fftw_execute(plan);
			// Do rescaling
			for (size_t i=0;i<N[0];i++) outp[i]/=static_cast<complex<double> >(N[0]);
		}
		else if (dim==2){
			fftw_plan plan=fftw_plan_dft_2d(N[0],N[1],reinterpret_cast<fftw_complex*>(&inp[0]),reinterpret_cast<fftw_complex*>(&outp[0]),FFTW_BACKWARD,FFTW_ESTIMATE);
			fftw_execute(plan);
			// Do rescaling
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++){
					outp[i*N[1]+j]/=static_cast<complex<double> >(N[0]*N[1]);
				}
		}
		else if (dim==3){
			fftw_plan plan=fftw_plan_dft_3d(N[0],N[1],N[2],reinterpret_cast<fftw_complex*>(&inp[0]),reinterpret_cast<fftw_complex*>(&outp[0]),FFTW_BACKWARD,FFTW_ESTIMATE);
			fftw_execute(plan);
			// Do rescaling
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++)
					for (size_t k=0;k<N[2];k++){
						outp[i*N[1]*N[2]+j*N[2]+k]/=static_cast<complex<double> >(N[0]*N[1]*N[2]);
				}
		}
	}
	else{
		cout << "Direction of FFTW is not defined." <<endl<<"1:Forward, 0:Backward"<<endl;abort();
	}

	fftw_cleanup_threads();
}

/*
 * dft_coeffs
 * It builds up coefs array.
 *
 *    if (dble(i0) <= n(1)/2.0d0) then
   coef = coef -(2.0d0*PI/Lx*dble(i0))**2
   else
   coef = coef -(2.0d0*PI/Lx*dble(i0-n(1)))**2
    endif
 */
void dft_coeffs(complex<double>* coefs,unsigned int order, size_t* N,double* L){

	double pi = 3.1415927;
	double	factors[3];	//	This one stores frequently used values.
	for (size_t i=0;i<3;i++){
		factors[i] = pow(2.0*pi/L[i],2);
	}

	if (order == 2){
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t k=0;k<N[2];k++){
					size_t indx = i*N[1]*N[2]+j*N[2]+k;
					coefs[indx] = 0;
					if (i<=(N[0]/2)){
						coefs[indx] -= static_cast<complex<double> >(factors[0] * pow(static_cast<double>(i),2.0));
					}
					else{
						coefs[indx] -= static_cast<complex<double> >(factors[0] * pow(static_cast<double>(i)-static_cast<double>(N[0]),2.0));
					}
					if (j<=(N[1]/2)){
						coefs[indx] -= static_cast<complex<double> >(factors[1] * pow(static_cast<double>(j),2.0));
					}
					else{
						coefs[indx] -= static_cast<complex<double> >(factors[1] * pow(static_cast<double>(j)-static_cast<double>(N[1]),2.0));
					}
					if (k<=(N[2]/2)){
						coefs[indx] -= static_cast<complex<double> >(factors[2] * pow(static_cast<double>(k),2.0));
					}
					else{
						coefs[indx] -= static_cast<complex<double> >(factors[2] * pow(static_cast<double>(k)-static_cast<double>(N[2]),2.0));
					}
				}
	}
	else{
		cout << "currently, order 1 is not supported. Aborting..."<<endl;
		abort();
	}


}


}


