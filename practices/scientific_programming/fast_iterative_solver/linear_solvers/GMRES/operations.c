
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mmio.h"
#include "operations.h"


double coo_to_csr(int* Ic, int* I, int* Jc, int* J_csr, int* J, double* valc, double* val, int M, int nz, int symm)
{
int n_csr=1, i, j, l, nz2=0;
if (symm == 0)
{
Ic[0]=I[0];
Jc[0]=J[0];
valc[0]=val[0];

// i: row number. j: nonzero element index
for (i=0; i< M; i++)
{
	for (j=1; j<nz; j++)
{
	if (I[j] == i)
{
Ic[n_csr] = i;
Jc[n_csr] = J[j];
valc[n_csr] = val[j];

n_csr++;
}
}
}

printf("Finished converting COO into COO2. \n");


    /* Converting COO2 format into CSR format */
printf("Starting converting COO2 into CSR. \n");
J_csr[0]=0;
J_csr[M] = nz;
n_csr=1;
    for (i=1; i<nz; i++)
    {
// we need to build up I_csr until it reaches to the last row
	if (Ic[i] > Ic[i-1])
{		J_csr[n_csr] = i;
		n_csr++;
// If we reached the last row, we can stop it
	if (Ic[i]== M)
	break;
}
   }

nz2 = nz;
}
// symmetric code starts here
else if (symm == 1)
{
J_csr[0] = 0;
n_csr = 0;
// searching for ith row & jth column
for (i=0; i<M; i++)
{
for (j=0; j<M; j++)
{
for (l=0; l< nz; l++)
{

if (((I[l] == i) && (J[l] == j)) || ((I[l] == j) && (J[l] == i)))
{
Jc[nz2] = j;
valc[nz2] = val[l];

if (Jc[nz2] <= Jc[nz2-1]) // bigger or equal to the previous value
{
J_csr[n_csr] = nz2;
n_csr++;
}
nz2++;
}
}
}
}
J_csr[n_csr] = nz2;

}// end symmetric code

return nz2;
}



void matrix_vector(int* J_csr, int* Jc, double* valc, double* Xv, double* y, int M)
{
//                CSR form, Jc: column, value vector, Xv: input vector, y: output, M: number of rows
int i, j;

// Initializing y where the result will be installed;

for (i=0; i<M; i++)
{
y[i] = 0;
}


for (i=0; i<M; i++)
	for (j=J_csr[i]; j<J_csr[i+1]; j++)
{
/*
printf("i, j %d %d \n", i, j);
printf("valc[j], Xv[Jc[j]] Jc[j] %e %e %d\n", valc[j], Xv[ Jc[j] ], Jc[j]);
*/
y[i] = y[i] + valc[ j ]*Xv[ Jc[j] ];

}
}



double norm(double* v, int k)
{
int i;
double s = 0;
for (i=0; i<k; i++)
{
s = s + v[i]*v[i];
}
s = sqrt(s);
return s;
}

double dot_p(double* a, double* b, int k)
{
//         a, b: input, c: output, k: size of input
int i;
double s = 0;
for (i=0; i<k; i++)
{
s = s+ a[i]*b[i];
}
return s;
}



double getkrylov(int* I, int* J, double* val, int* Jc, int* J_csr, double* valc, double* H, int* H_csc, double* V, double* V_c,  double* w, const int j, int M, int N, int nz, double* g, double* c, double* s, int symm, int* precond, int* precond_csr, double* precond_val, int precondition)//  nz could be needed
{

int i, k, l;
double *temp_arr = (double *) calloc(M, sizeof(double)); // temporary array
for (i=0; i<M; i++)
{
V_c[i] = V[M*i + j];
// I initialized w[i] and it works.
w[i]=0;
}

// calculate w = A * V(:,j)
if (precondition == 0)
matrix_vector(J_csr, Jc, valc, V_c, w, M);
else if (precondition == 1)
{
matrix_vector(J_csr, Jc, valc, V_c, temp_arr, M);
matrix_vector(precond_csr, precond, precond_val, temp_arr, w, M);
}
else if (precondition == 2)
{
	int precond_count = M*(M+1)/2;
matrix_vector(J_csr, Jc, valc, V_c, temp_arr, M);
H_inv_g(precond_val, precond_csr, M, precond_count, temp_arr, w);
}

int index = H_csc[j]; // starting point of jth column.
double temp1, temp2, temp3;

// working on ith rows
for (i=0; i<j+1; i++)
{
// copy i-th column
for (k=0; k<M; k++)
{
// V_c is copy of i-th column
V_c[k] = V[M*k + i];
}

// H ith row in jth column
H[index+ i] = dot_p(V_c, w, M);


for (l=0; l<M; l++)
{
w[l] = w[l] - H[index + i]*V_c[l];
}
}
H[index+j+1] = norm(w, M);

// storing next basis
if (j < M-1)
{
for (i=0; i<M; i++)
{
V[i*N + j + 1] = w[i]/H[index+j+1];
//printf("next basis V[%d] = %e\n\n",i*N+j+1, V[i*N + j + 1]);
}
}

for (k=1; k <= j; k++)
{
index = H_csc[j];
temp1= c[k-1]*H[index + k-1] + s[k-1]*H[index + k];
temp2= -s[k-1]*H[index + k-1] + c[k-1]*H[index + k];
H[index + k-1] = temp1;
H[index + k] = temp2;
}

temp3= sqrt( (H[index + j])*(H[index + j]) + (H[index + j+1])*(H[index + j+1]) );
c[j] = H[index + j] / temp3;
s[j] = H[index + j + 1] / temp3;


temp1 = c[j]*H[H_csc[j]+j] + s[j]*H[H_csc[j]+j+1];
temp2 = -s[j]*H[H_csc[j]+j] + c[j]*H[H_csc[j]+j+1];


H[index+j] = temp1;
H[index+j+1] = temp2;

g[j+1] = -s[j]*g[j];
g[j] = c[j]*g[j];
//printf("g[%d] = %e\n",j+1,g[j+1]);


free(temp_arr);

return g[j+1];
}



void H_inv_g(double* H, int* H_csc, int M, int H_size, double* g, double* h_inv_g)
{//          Hessian mat,   CSC form, #row, size of H,   g vector,  output
int i, j;
double s=0;

// printf("\n\n starting H_inv_g \n\n");

/*
// output initialization
for (i=0; i<M; i++)
{
h_inv_g[i] = 0;
}
*/
h_inv_g[M-1] = 1/H[H_csc[M]-1]*g[M-1];

for (i =M-2; i >= 0; i--)
{
s= 0;
for (j = i+1; j <= M-1; j++)
{
//index = H_csc[j] + i;
s= s+ H[H_csc[j] + i]*h_inv_g[j];
/*printf("index = %d\n", index);
printf("h_inv_g[%d] = %e\n", j, h_inv_g[j]);*/
}
//printf("H[i] + i = %d");
h_inv_g[i] = (g[i] - s)/(H[H_csc[i]+i]);
}
}

/*
correct ans
   2.195133536330757
   0.360301879288834
  -0.227093182321466
*/




void dense_matrix_vector(double* A, int M, int N, double* x, double* y)
{ //                    inp Matrix    row   column   inp      output
int i, j;

// Initializing y where the result will be installed;

for (i=0; i<M; i++)
{
y[i] = 0;
}

// running program
for (i=0;i<M; i++)
for (j=0; j<N; j++)
y[i] = y[i] + A[N*i+ j]*x[j];
/*
printf("\nM = %d \n",M);
// output
printf("\ny\n");
for (j=0; j<M; j++)
printf("%e\n", y[j]);
*/
}


int Precondition(int* Jc, int* J_csr, double* valc, int M, int nz, int* precond, int* precond_csr, double* precond_val, int precondition)
{           // 1: input: CSR format,  value array, rank(A), nz,  output: column, csr,            value,           sizeof array
			// 2: input: original I, J,    val,  rank(A), nz,  output: null,  CSC,           value,          sizeof array
	int i, j;
int precond_count = 0;
int *I = Jc;
int *J = J_csr;
double *val = valc;

printf("Jacobian preconditioner will be applied\n");

if (precondition == 1)  // Jacobi
{
for (i=0; i<M; i++)
{
for (j=J_csr[i]; j<J_csr[i+1]; j++)
{
if (Jc[j] == i)
{
precond[i] = i;
precond_val[i] = 1/valc[j];
precond_csr[i] = i;
precond_count++;
}
}
}
precond_csr[M] = precond_count;
}


else if (precondition == 2) // Gauss-Seidel method. I am making upper-triangular matrix in CSC format.
{	// 2: input: original I, J,    val,  rank(A), output: null,  CSC,           value,          sizeof array

int sizeof_precond;
//	FILE * out2;
//	int count1 = 0, count2 = 0;

	printf("Gauss-Seidel preconditioner will be applied\n");


	// initialization of precond_csr
for (i=0; i<=M; i++)
{
	precond_csr[i] = i*(i+1)/2;
}
sizeof_precond = precond_csr[M];

// Saving nonzero terms into precond_val
for (i=0; i< nz; i++)
{
if (I[i] <= J[i])
{
precond_val[precond_csr[J[i]] + I[i]] = val[i];

if ((I[i] == J[i]) && (val[i] == 0))
{
	printf("Zero diagonal term");
	exit(-1);
}

}
//printf("precond_count = %d \n",precond_count);


/*   This calloc below yield a serious memory problem*/

// Still need to implement inverse of precond
// It maybe necessary to calculate inverse of (L+D) when precondition == 1;
}
precond_count = sizeof_precond;
}

return precond_count;
}

