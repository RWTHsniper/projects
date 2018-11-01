

/*
 * myfun.h
 *
 *  Created on: Jul 15, 2016
 *      Author: jaeyong
	It has my library of math functions
 */

#ifndef MYFUN_H_
#define MYFUN_H_


void Inv_sign(double*A, int M);
void swapbytes(char* array, int nelem, int elemsize);
void address_swaper(int **A, int **B);
void address_swaper(double **A, double **B);
void read_data(double* data, size_t nn, unsigned entries, char* file_name);
void write_data(double* data, size_t nn, unsigned entries, char* file_name);
void Init(double*A, int M, double N); 	// Initialize arrays, size of array, initializing value
double dot_p(double *A, double *B, int M);
void M_transpose(double *B, double *A, int M);
double determinant(double *A, int M);
void check_value(double *A, int M); // checks every terms in an array
void D_matrix_vector(double* y, double* A, int M, int N, double* x);
double Inv_M(double* result, double* A, int M);
void D_print(double* A, int M, int N);
void I_print(int* A, int M, int N);
void M_fprint(char* name, double* A, int M, int N); // print M by N matrix on text file
void M_I_fprint(char* name, int* A, int M, int N); // print M by N matrix on text file
double abs_max(double *A, int M); // find abs_max
double norm(double* A, int M);
void Unit_v(double* A, double* B, int M);
double M_inf_norm(double* A, int r, int c);
void scale(double* A, double* B, int M, double factor); // store A =B*factor
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

// CSR section
double coo_to_csr(int* Ic, int* I, int* Jc, int* J_csr, int* J, double* valc, int *I_csc, double* val, int M, int nz, int symm);
void CSR_matrix_vector(int* J_csr, int* Jc, int* I_csc, int* I, double* val, double* valc, double* Xv, int symm, double* y, int M);
void CSR_vector_matrix(int* J_csr, int* Jc, int* I_csc, int* I, double* val, double* valc, double* Xv, int symm, double* y, int M);
void CSR_CG_solver(double *x, int *J_csr, int *Jc, int *I_csc, int *I, double *val, double *valc, int symm, double *Y, int M, double tol);
double power_iteration(int *I_csc, int *I, int *J_csr, int *Jc, double *val, double *valc, int symm, int M, int N, double tol);
double inv_power_iteration(int *I_csc, int *I, int *J_csr, int *Jc, double *val, double *valc, int symm, int M, int N, double tol);
double CSR_lanczos(int *I_csc, int *I, int *J_csr, int *Jc, double *val, double *valc, int symm, int M, int N);


#endif /* MYFUN_H_ */




