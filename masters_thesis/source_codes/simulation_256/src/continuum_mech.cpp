/*
 * continuum_mech.cpp
 *
 *  Created on: Nov 1, 2017
 *      Author: jaeyong
 */

#include "continuum_mech.h"

using namespace std;

namespace finite_strain{

}

namespace elasticity{

void convert_P_2_Sig(double* sig, double*P,double* F){
//	sig = 1.0/det(F)*P'*F;
	double detF = math_op::determinant(F,3);
	for (size_t i=0;i<3;i++)
		for (size_t j=0;j<3;j++){
			sig[i*3+j]=0.0e0;
			for (size_t k=0;k<3;k++){
				sig[i*3+j]+= P[k*3+i]*F[k*3+j];
			}
			sig[i*3+j]/=detF;
		}
}

namespace hyper{
namespace neo_hook{
void piola_kirch_stress(double* S, double* F,double mu, double lam,unsigned int order){

	double C[3][3], invC[3][3];
	finite_strain::right_cauchy_green_T(&C[0][0],F);
	double detC = math_op::determinant(&C[0][0],3);
	math_op::Inv_M(&invC[0][0],&C[0][0],3);

// 2nd order Piola-Kirchhoff
		for (size_t i=0;i<3;i++)
			for (size_t j=0;j<3;j++){
				double dij = double(i==j);	// Kroneker's delta
				S[i*3+j]=mu*(dij-invC[i][j])+lam/2.0*(detC-1.0)*invC[i][j];
//				S[i*3+j]=mu*(eye[i][j]-invC[i][j])+lam/2.0*(detC-1.0)*invC[i][j];
			}
	if (order == 1){// First order
//	P = F*S
		double P[3][3];
		for (size_t i=0;i<3;i++)
			for (size_t j=0;j<3;j++){
				double tmp = 0;
				for (size_t k=0;k<3;k++){
					tmp+=F[i*3+k]*S[k*3+j];
				}
//				cout << tmp<<endl;
				P[i][j]=tmp;
			}
//		Load calculated values onto output array
		for (size_t i=0;i<3;i++)
			for (size_t j=0;j<3;j++)	S[i*3+j]=P[i][j];
	}
	else if (order !=2){
		cout << "order of Piola-Kirchhoff tensor should be defined" << endl;
		abort();
	}
	}


}
}


}
