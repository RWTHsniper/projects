/*
 * helpers.h
 *
 *  Created on: Oct 27, 2017
 *      Author: jaeyong
 */

#ifndef HELPERS_H_
#define HELPERS_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

namespace helpers{
//	This is for multidimensional array
template<typename T>
void print_matrix(T* A, unsigned int m, unsigned int n){

	for (size_t i=0; i<m; i++){
		for (size_t j=0; j<n; j++){
			std::cout << A[i*n+j] <<" ";
		}	std::cout << std::endl;
	}

}
//	This is for multidimensional vector
template<typename T>
void print_matrix(T& A){
	if (A.size()==0)	cout << "Size of vector is not defined!!!. Impossible to print" << endl;
		for (size_t i=0;i<A.size();i++){
			for (size_t j=0;j<A[i].size();j++){
				cout << A[i][j] << " ";
			}			cout << endl;
		}
}

//	I want to make multi-dimensional vector creator.
//	2D vector
template<typename T>
void create_multidim_vector(vector<vector<T> >& A,size_t m,size_t n, T initi){

	A.resize(m);
for (size_t i=0;i<m;i++){
	A[i].resize(n);
	for (size_t j=0;j<n;j++){
		A[i][j] = initi;
	}
}
}

// 3D vector
template<typename T>
void create_multidim_vector(vector<vector<vector<T> > >& A,size_t m,size_t n,size_t o, T initi){

	A.resize(m);
for (size_t i=0;i<m;i++){
	A[i].resize(n);
	for (size_t j=0;j<n;j++){
		A[i][j].resize(o);
		for (size_t k=0;k<o;k++){
			A[i][j][k] = initi;
		}
	}
}
}
// 4D vector
template<typename T>
void create_multidim_vector(vector<vector<vector<vector<T> > > >& A,size_t m,size_t n,size_t o,size_t p,T initi){

	A.resize(m);
for (size_t i=0;i<m;i++){
	A[i].resize(n);
	for (size_t j=0;j<n;j++){
		A[i][j].resize(o);
		for (size_t k=0;k<o;k++){
			A[i][j][k].resize(p);
			for (size_t l=0;l<p;l++){
				A[i][j][k][l] = initi;
				}
			}
		}
	}
}

string put_leading_zeros(unsigned int width,unsigned int number);

}

#endif /* HELPERS_H_ */
