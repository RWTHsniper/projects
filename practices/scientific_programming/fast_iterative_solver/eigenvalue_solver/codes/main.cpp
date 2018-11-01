/*
   ============================================================================
Name        : fis_p3_1.c
Author      : jaeyong
Version     :
Copyright   : Your copyright notice
Description : Hello World in C, Ansi-style
============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "myfun.h"
#include "mmio.h"
#include <time.h>

#include <iostream>
#include <sstream>
#include <string>


using namespace std;


int main(int argc, char *argv[]) {

	int ret_code;
    FILE *f;
    MM_typecode matcode;
    int M, N, nz, symm = 0;


    // Input read and checking procedure.

    if (argc < 2)
    {
	fprintf(stderr, "Usage: %s [matrix-market-filename]\n", argv[0]);
	exit(1);
    }
    else
    {
	if ((f = fopen(argv[1], "r")) == NULL)
	    exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
	printf("Could not process Matrix Market banner.\n");
	exit(1);
    }

    //  This is how one can screen matrix types if their application 
    //  only supports a subset of the Matrix Market data types.      

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
	    mm_is_sparse(matcode) )
    {
	printf("Sorry, this application does not support ");
	printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
	exit(1);
    }

    // find out size of sparse matrix 

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
	exit(1);


    printf("Input check is finished. \n");




    // basic arrays for CSC/CSR.
    int *I = (int *) calloc(nz, sizeof(int));
    int *Ic = (int *) calloc(nz, sizeof(int));
    int *J_csr = (int *) calloc((M+1), sizeof(int));
    int *I_csc = (int *) calloc((M+1), sizeof(int));
    int *J = (int *) calloc(nz, sizeof(int));
    //	int *Jc = (int *) calloc(2*nz, sizeof(int));
    int *Jc = (int *) calloc(nz, sizeof(int));
    double *val = (double *) calloc(nz, sizeof(double));
    //	double *valc = (double *) calloc(2*nz, sizeof(double));
    double *valc = (double *) calloc(nz, sizeof(double));
    double *y = (double *) calloc(N, sizeof(double));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    // I: X, J: Y reading data
    for (int i=0; i<nz; i++)
    {
	fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
	I[i]--;  /* adjust from 1-based to 0-based */
	J[i]--;
    }

    if (f !=stdin) fclose(f);

    // remind that the initial data are arranged in columns
    // I will compress J


    if (mm_is_symmetric(matcode) == 1)
    {
	symm = 1;
    }

    nz = coo_to_csr(Ic, I, Jc, J_csr, J, valc, I_csc, val, M, nz, symm);

	// finishing reading CSR format matrix



double time;
double lambda;
double t_lambda = 9.5986080894852857e+03;

	time =clock();
    lambda = power_iteration( I_csc, I, J_csr, Jc, val, valc, symm, M, N, 1e-6);
	printf("Power iteration time = %e \n", (clock()-time)/CLOCKS_PER_SEC);
	printf("s3rmt3m3 error %e\n", fabs(t_lambda-lambda));

    //inv_power_iteration( I_csc, I, J_csr, Jc, val, valc, symm, M, N, 1e-8);

	time =clock();
    lambda = CSR_lanczos(I_csc, I, J_csr, Jc, val, valc, symm, M, N);
	printf("Lanczos time = %e \n", (clock()-time)/CLOCKS_PER_SEC);
	printf("s3rmt3m3 error %e\n", fabs(t_lambda-lambda));


    free(I);
    free(Ic);
    free(J_csr);
    free(I_csc);
    free(J);
    free(Jc);
    free(val);
    free(valc);
    free(y);



    return EXIT_SUCCESS;
}
