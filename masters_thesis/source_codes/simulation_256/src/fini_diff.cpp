/*
 * fini_diff.cpp
 *
 *  Created on: Nov 10, 2017
 *      Author: jaeyong
 */

#include "fini_diff.h"

using namespace std;

typedef Eigen::Triplet<double> T;

namespace FD{

size_t	global_indx(long int i,long int j,long int k,size_t* N){
	size_t indx = 0;
	long int lN[3];	//	long N
	for (size_t i=0;i<3;i++)	lN[i] = static_cast<long int>(N[i]);
	if (i<0){
		i += lN[0];
	}
	else if (i>=lN[0]){
		i -= lN[0];
	}
	if (j<0){
		j += lN[1];
	}
	else if (j>=lN[1]){
		j -= lN[1];
	}
	if (k<0){
		k += lN[2];
	}
	else if (k>=lN[2]){
		k -= lN[2];
	}
	indx = static_cast<size_t>(i*lN[1]*lN[2]+j*lN[2]+k);
	return indx;
}

/*
! x direction
        call global_index_mapper(g_indx, i-1,j,k,n )
	indx_2 = g_indx
        dA11(indx_2) = dA11(indx_2) -con(1)
        call global_index_mapper(g_indx, i+1,j,k,n )
 	indx_2 = g_indx
        dA11(indx_2) = dA11(indx_2) -con(1)
 ! y direction
        call global_index_mapper(g_indx, i,j-1,k,n )
	indx_2 = g_indx
        dA11(indx_2) = dA11(indx_2) -con(2)
        call global_index_mapper(g_indx, i,j+1,k,n )
	indx_2 = g_indx
        dA11(indx_2) = dA11(indx_2) -con(2)
 ! z direction
        call global_index_mapper(g_indx, i,j,k-1,n )
	indx_2 = g_indx
        dA11(indx_2) = dA11(indx_2) -con(3)
        call global_index_mapper(g_indx, i,j,k+1,n )
	indx_2 = g_indx
        dA11(indx_2) = dA11(indx_2) -con(3)
!print*,'dA11 check'
!print*,dA11(1:u_l)

!        print*,'Lets call CSR build up'
call load_vector_CSR(row_indx,col,values,dA11,indx_ijk,int_nxnynz,int_nxnynz,sizeofvalues,1)
 *
 */
/*
void	build_CD_matrix(spar::spar_mat& fdmat,size_t* N,double* dx,double alph, unsigned int model){
	cout << "dx "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<endl;
	cout << "N "<<N[0]<<" "<<N[1]<<" "<<N[2]<<endl;
	size_t inxnynz = N[0]*N[1]*N[2];
	vector<double> con(3,0);
	for (size_t i=0;i<3;i++){
		if (N[i]==1){
			con[i] = 0.0;
		}
		else{
			con[i] = alph/(dx[i]*dx[i]);
		}
	}
	if (model==0){

		for (size_t i=0;i<N[0];i++){
			long int ii= static_cast<long int>(i);
			for (size_t j=0;j<N[1];j++){
				long int jj=static_cast<long int>(j);
				for (size_t k=0;k<N[2];k++){
					size_t indx_ijk =i*N[1]*N[2]+j*N[2]+k;
					vector<double>	rowvect(inxnynz,0.0);
					rowvect[indx_ijk] = 1.0e0+2.0e0*(con[0]+con[1]+con[2]);
					long int kk = static_cast<long int>(k);
					//	x direction
					size_t gindx = global_indx(ii-1,jj,kk,N);
					rowvect[gindx] -= con[0];
					gindx = global_indx(ii+1,jj,kk,N);
					rowvect[gindx] -= con[0];
					//	y direction
					gindx = global_indx(ii,jj-1,kk,N);
					rowvect[gindx] -= con[1];
					gindx = global_indx(ii,jj+1,kk,N);
					rowvect[gindx] -= con[1];
					//	z direction
					gindx = global_indx(ii,jj,kk-1,N);
					rowvect[gindx] -= con[2];
					gindx = global_indx(ii,jj,kk+1,N);
					rowvect[gindx] -= con[2];
					fdmat.build_csr(&rowvect[0]);
				}
			}
		}
	}
	else{
		cout << "Other models like higher order CD is not supported yet."<<endl;
		abort();
	}
}
*/
/*	build matrix of Eigen
 *
 */
void	build_CD_matrix(Eigen::SparseMatrix<double,Eigen::RowMajor>& argMat,size_t* N,double* dx,double alph, unsigned int model){
	cout << "dx "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<endl;
	cout << "N "<<N[0]<<" "<<N[1]<<" "<<N[2]<<endl;
	size_t inxnynz = N[0]*N[1]*N[2];
	vector<double> con(3,0);
	size_t dim=0;
	for (size_t i=0;i<3;i++){
		if (N[i]==1){
			con[i] = 0.0;
		}
		else{
			con[i] = alph/(dx[i]*dx[i]);
			dim++;
		}
	}
	if (model==0){
		std::vector<Eigen::Triplet<double> > tripletlist;
		tripletlist.reserve(inxnynz*(1+2*dim));

		for (size_t i=0;i<N[0];i++){
			long int ii= static_cast<long int>(i);
			for (size_t j=0;j<N[1];j++){
				long int jj=static_cast<long int>(j);
				for (size_t k=0;k<N[2];k++){
					size_t indx_ijk =i*N[1]*N[2]+j*N[2]+k;
					vector<double>	rowvect(inxnynz,0.0);
					rowvect[indx_ijk] = 1.0e0+2.0e0*(con[0]+con[1]+con[2]);
					long int kk = static_cast<long int>(k);
					//	x direction
					size_t gindx = global_indx(ii-1,jj,kk,N);
					rowvect[gindx] -= con[0];
					gindx = global_indx(ii+1,jj,kk,N);
					rowvect[gindx] -= con[0];
					//	y direction
					gindx = global_indx(ii,jj-1,kk,N);
					rowvect[gindx] -= con[1];
					gindx = global_indx(ii,jj+1,kk,N);
					rowvect[gindx] -= con[1];
					//	z direction
					gindx = global_indx(ii,jj,kk-1,N);
					rowvect[gindx] -= con[2];
					gindx = global_indx(ii,jj,kk+1,N);
					rowvect[gindx] -= con[2];
//					fdmat.build_csr(&rowvect[0]);
					for (size_t iter=0;iter<N[0]*N[1]*N[2];iter++){
					if (fabs(rowvect[iter])>0.0)	{
//						argMat.insert(iter,indx_ijk)= rowvect[iter];
						tripletlist.push_back(T(iter,indx_ijk,rowvect[iter]));
					}
					}
				}
			}
		}
		argMat.setFromTriplets(tripletlist.begin(), tripletlist.end());
	}
	else{
		cout << "Other models like higher order CD is not supported yet."<<endl;
		abort();
	}
}


}
