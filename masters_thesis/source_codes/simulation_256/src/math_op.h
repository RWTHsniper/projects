
#ifndef MATH_OP_H_
#define MATH_OP_H_


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <omp.h>
#include "helpers.h"
//#include "/usr/include/fftw3.h"

namespace math_op{

template <typename T>
inline void initiali_arr(T *A,T val, size_t M)
{
    for (size_t i=0; i<M; i++)
    {
	A[i] = val;
    }
}

template <typename T>
inline void copy_arr(T *A,T* B, size_t M)
{
    for (size_t i=0; i<M; i++)
    {
	A[i] = B[i];
    }
}

template <typename T>
inline void exchange(T *A, T *B, size_t M )
{
    T	C;
    for (size_t i=0; i<M; i++)
    {
	C = A[i];
	A[i] = B[i];
	B[i] = C;
    }
}

template <typename T>
T determinant(T* A, size_t M)
{
	size_t n = M;
    T tol = 1e-20;
    T det=1;
    T (*B) = (T(*)) calloc(M*M, sizeof(T));
    if(n==1) return A[0*M + 0];
    if(n==2) return (A[0*M + 0]*A[1*M + 1] - A[0*M + 1]*A[1*M + 0]);
    if(n==3) return (A[0*M + 0]*((A[1*M + 1]*A[2*M+2]) - (A[2*M+1]*A[1*M+2])) -A[0*M+1]*(A[1*M+0]*A[2*M+2] - A[2*M+0]*A[1*M+2]) + A[0*M+2]*(A[1*M+0]*A[2*M+1] - A[2*M+0]*A[1*M+1]));


    // Do row reduction to make echelon form
    for (size_t i=0;i<M*M;i++) B[i]=A[i];


    //printf("trying to make echelon form\n");
    for (size_t k=0; k<M-1; k++) // k: reference row
    {
	for (size_t i=k+1; i<M; i++) // i: target row
	{

	    if (fabs(B[k*M+k]) < tol)
	    {
		//printf("[%d %d]exchaning rows\n",i,k);
		for (size_t l=k+1; l<M; l++)
		{
		    if (fabs(B[l*M+k]) > tol)
		    {
			det = det*-1;
			exchange(&(B[k*M+k]), &(B[l*M+k]), (M-k) );
			break;
		    }
		}
	    }

	    if (fabs(B[k*M+k]) < tol)
	    {

		std::cout << "singular pivot ,"<<i<<", "<<k<<", "<<M<<", "<<B[i*M+k]<<", "<<B[k*M+k]<<", "<< std::endl;
		printf("B printing \n\n");
		helpers::print_matrix(B, M, M);

		return 0;
	    }

	    T fact = B[i*M+k]/B[k*M+k];
	    //printf("[%d, %d] %e\n",i, k, fact);
	    for (size_t j=k; j<M; j++) // calculating jth column
	    {
		B[i*M+j] = B[i*M+j] -fact*B[k*M+j];
	    }


	}
    }

    for (size_t i=0; i<M; i++)
    {
	det = det*B[i*M+i];
    }

    if(isfinite(det) == 0)
    {
	printf("\nDeterminant is fucked up as you can see %e \n",det);
	exit(-1);
    }

    free (B);
    return det;
}

template <typename T>
inline void TP(T* out, T* inp, size_t m, size_t n){
/* *
 * inp[m][n]
 * out[n][m]
 * */
	if (&out[0]==&inp[0]){
//		printf("same inp & outp array in TP\n");
	    T * tmp = (T(*)) calloc((n)*(m), sizeof(T));
		for (size_t i=0;i<m;i++)
			for(size_t j=0;j<n;j++){
				// tmp[j][i] = inp[i][j]
				tmp[j*m+i]=inp[i*n+j];
			}
		for (size_t i=0;i<n;i++)
			for(size_t j=0;j<m;j++){
				// out[i][j]=tmp[i][j]
				out[i*m+j]=tmp[i*m+j];
			}
	    free (tmp);
	}
	else{
		for (size_t i=0;i<m;i++)
			for(size_t j=0;j<n;j++){
				// out[j][i] = inp[i][j]
				out[j*m+i]=inp[i*n+j];
			}
	}

}

template <typename T>
T Inv_M(T *result, T *A,size_t M) // result storage, input matrix, size of matrix
    // input and output should be different!!!!
{
//    printf("Inversion procedure start!\n");

    // M is dimension of input matrix A

    T *C = (T(*)) calloc((M-1)*(M-1), sizeof(T));
    T tol = 1e-12;
    T det=0;


    det = determinant(A, M);

    if (fabs(det) < tol)
    {
	printf("\nsingular matrix. Fucked up.\n");
	printf("tolerance = %e\n",tol);
	helpers::print_matrix(A,M,M);
	exit(-1);
    }

    // Build matrix of minors
    for (size_t i=0; i<M; i++)
    {

	for (size_t j=0; j<M; j++)
	{
		size_t count=0;
	    // cofactor amtrix for B[i][j]
	    for (size_t k=0; k<M; k++)
	    {
		for (size_t l=0; l<M; l++)
		{
		    if ((k == i) || (l == j))
		    {
			continue;
		    }
		    else
		    {
			C[count] = A[k*M + l];
			count++;
		    }
		}
	    }
	    result[i*M + j] = determinant(C,M-1);
	    result[i*M + j] = result[i*M + j]/det;
	    if (((i+j) % 2) == 1)
	    {
		result[i*M + j] = -result[i*M + j];
	    }
	}
    }

    math_op::TP(result, result, M,M);


    free (C);


    T inv_check= determinant(result,M);
    if ((1 - tol >det*inv_check) || (1 + tol < det*inv_check))
    {
	printf("\n[%e %e]Inverse matrix calculation error\n", det, inv_check);
	printf("tolerance = %e\n",tol);
	printf("result printing \n\n");
	helpers::print_matrix((result), M, M);
	abort();
    }

    return det;
}

template <typename T>
inline T trace(T*A, size_t m){
T outp = 0;
for (size_t i=0;i<m;i++){
	outp += A[i*m+i];
}
return outp;
}

template <typename T>
inline void add(T*C, T* A, T* B,T coef1,T coef2, size_t len){
//	C = A + B
	for (size_t i=0;i<len;i++){
		C[i]=coef1*A[i]+coef2*B[i];
	}
}

template <typename T>
inline void scale(T* B,T* A,T coef,unsigned int len){//	B = coef*A;

	for (size_t i=0;i<len;i++){
		B[i]=coef*A[i];
	}

}

/*
 * dotp
 * calculate dot product between two array
 */
template <typename T>
inline T dotp(T* A,T* B,unsigned int n){//	B = coef*A;

	double outp = 0.0;
	for (size_t i=0;i<n;i++){
		outp += A[i]*B[i];
	}
	return outp;
}

/*
 * norm2:
 * calculate 2-norm of an array
 */
template <typename T>
inline T norm2(T* A,unsigned int n){//	B = coef*A;

	double outp = 0.0;
	for (size_t i=0;i<n;i++){
		outp += A[i]*A[i];
	}
	outp = sqrt(outp);
	return outp;
}

/*	eig_power
 * calculate largest absolute value of eigenvalue using power iteration.
 */
template <typename T>
inline T eig_power(T* A,unsigned int n){//	B = coef*A;

	vector<T>	q(n,0); q[0]=1.0;
	vector<T>	z(n,0);

	size_t maxiter = 2*n;
	T lam = 0;
	T lam_pr=0;
	T normz = 0;

	for (size_t it=0; it<maxiter;it++){
		//	calculate z = A*q
		for (size_t i=0;i<n;i++){
			z[i]=0;
			for (size_t j=0;j<n;j++){
					z[i]+= A[i*n+j]*q[j];
			}
		}
		normz = math_op::norm2(&z[0],n); lam = 0;
		for (size_t i=0;i<n;i++) q[i] = z[i]/normz;

		for (size_t i=0;i<n;i++)
			for (size_t j=0;j<n;j++){
				lam+=q[i]*A[i*n+j]*q[j];
			}
		if (lam/lam_pr -1 < 1.0e-5){
			break;
		}
		else	{
//			cout << "current estimation of lam = "<<lam<<endl;
			lam_pr = lam;
		}
	}

	return lam;
}


template <typename T>
//	Add 2D vectors
inline void add(vector<vector<T> >&C, vector<vector<T> >&A, vector<vector<T> >&B){

	if (C.size()==A.size()==B.size()) cout << "Dimension of vector arrays mismatch" << endl;

	for (size_t i=0;i<A.size();i++){
		for (size_t j=0;j<A[i].size();j++){
			C[i][j]=A[i][j]+B[i][j];
		}
	}
}


}



#endif
