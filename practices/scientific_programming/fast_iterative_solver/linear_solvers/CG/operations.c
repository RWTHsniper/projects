
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mmio.h"
#include "operations.h"


double coo_to_csr(int* Ic, int* I, int* Jc, int* J_csr, int* J, double* valc, int *I_csc, double* val, int M, int nz, int symm)
{
	int n_csr=1, i, j, nz2=0;
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
		printf("CSC, CSR operation \n");
		int count_csc=1;
		// build up I_csc
		I_csc[0] = 0;
		I_csc[M] = nz;
		for (i=1; i<nz; i++)
		{
			if (J[i] > J[i-1])
			{
				I_csc[count_csc] = i;
				count_csc++;
			}
		}



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
	}// end symmetric code

	return nz2;
}



void matrix_vector(int* J_csr, int* Jc, int* I_csc, int* I, double* val, double* valc, double* Xv, int symm, double* y, int M)
{
	//                CSR form, Jc: column, value vector, Xv: input vector, y: output, M: number of rows
	int i, j;

	// Initializing y where the result will be installed;
	for (i=0; i<M; i++)
	{
		y[i] = 0;
	}

	if (symm ==0)
	{
		for (i=0; i<M; i++)
			for (j=J_csr[i]; j<J_csr[i+1]; j++)
			{

				y[i] = y[i] + valc[ j ]*Xv[ Jc[j] ];
			}
	}
	else
	{
		// do calculation in row using csr
		for (i=0; i<M; i++)
			for (j=J_csr[i]; j<J_csr[i+1]; j++)
			{
				y[i] = y[i] + valc[ j ]*Xv[ Jc[j] ];
			}
/*
		printf("row operation\n");
		for (i=0;i<M;i++)
		{
			printf("Y[%d] = %e\n",i, y[i]); 
		}
*/
		// do calculation in column using csc

		for (i=0; i<M; i++)
		{
			for (j=I_csc[i] + 1; j<I_csc[i+1]; j++)
			{
				y[i] = y[i] + val[ j ]*Xv[I[j]];
			}
		}

/*		printf("column operation\n");
		for (i=0;i<M;i++)
		{
			printf("Y[%d] = %e\n",i, y[i]); 
		}
*/
	}// end else
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

	printf("\nM = %d \n",M);
	// output
	printf("\ny\n");
	for (j=0; j<M; j++)
		printf("%e\n", y[j]);

}

void desne_inverse_matrix_vector(double* A, int M, double* x, double*y)
{
	// It is presumed that A is upper triangular matrix
	int i, j, s;

	y[M-1] = 0;

	// i: row, j: column
	for (i=M-1; i>=0;i--)
	{
		s=0;
		for (j=i+1;j<M;j++)
		{
			s = s + A[i*M+j]*y[j];
		}
		y[i] = (x[i]-s)/A[i*M+i];
	}
}
