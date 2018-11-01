/*
 * continuum_mech.hpp
 *
 *  Created on: Nov 1, 2017
 *      Author: jaeyong
 */

#ifndef CONTINUUM_MECH_H_
#define CONTINUUM_MECH_H_

#include <stdio.h>
#include <stdlib.h>
#include "math_op.h"
#include "helpers.h"


namespace finite_strain{
inline void right_cauchy_green_T(double* C,double* F){

	for (size_t i=0;i<3;i++)
		for (size_t j=0;j<3;j++){
			C[i*3+j]=0.0;
			for (size_t k=0;k<3;k++){
			//	C = F'*F;
			C[i*3+j]+=F[k*3+i]*F[k*3+j];
			}
		}
}
inline void green_lagrangian_T(double* E,double* F){
	for (size_t i=0;i<3;i++)
		for (size_t j=0;j<3;j++){
			E[i*3+j]=0.0;
			for (size_t k=0;k<3;k++){
				//	E = 0.5*(C-I);
				E[i*3+j]+=(F[k*3+i]*F[k*3+j]);
			}
			E[i*3+j]-=double(i==j);
			E[i*3+j]*=0.5;
		}
}
}

namespace elasticity{
void convert_P_2_Sig(double* sig, double*P,double* F);

namespace hyper{
namespace neo_hook{
void piola_kirch_stress(double* S, double* F,double mu, double lam,unsigned int order);
}
}

}




#endif /* CONTINUUM_MECH_H_ */
