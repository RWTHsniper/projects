/*
 * myfun.h
 *
 *  Created on: Jul 15, 2016
 *      Author: jaeyong
	It has my library of math functions
 */

#ifndef MYFUN_H_
#define MYFUN_H_


void Initialization(double* u_s, double* f, int N);
void swapbytes(char* array, int nelem, int elemsize);
void read_data(double* data, size_t nn, unsigned entries, char* file_name);
void write_data(double* data, size_t nn, unsigned entries, char* file_name);
void Inv_sign(double*A, int M);
void Init(double*A, int M, double N); 	// Initialize arrays, size of array, initializing value
void M_transpose(double *B, double *A, int M);
double determinant(double *A, int M);
void check_value(double *A, int M); // checks every terms in an array
void D_matrix_vector(double* y, double* A, int M, int N, double* x);
double Inv_M(double* result, double* A, int M);
void M_print(double* A, int M, int N);
void M_fprint(char* name, double* A, int M, int N); // print M by N matrix on text file
double abs_max(double *A, int M); // find abs_max
double norm(double* A, int M);
double M_inf_norm(double* A, int r, int c);
void scale(double* A, int M, double factor);
void copy(double* B, double* A, int M); // Copy A onto B
void exchange(double *A, double *B, int M );
void add(double* C, double* A, double* B, int M);
void subtract(double* C, double* A, double* B, int M); // A -B = C, M:length of array
void GS(double* u0, double* f, int nu, int N);
void Restriction(double* u_c, double* u, int N);
void Prolongation(double* u, double* u_c, int N);
void Laplacian(double* u_t, double* u, double h, int N);
double r_inf_norm(double *u, double* f, int N);
void Inv_Laplacian(double* result, double* u, double h, int N);
void MG(int l, double* u, double* f, int N, int gamma, int nu1, int nu2);


#endif /* MYFUN_H_ */
