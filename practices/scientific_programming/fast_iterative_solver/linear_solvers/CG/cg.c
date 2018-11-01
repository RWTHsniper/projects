/*
 *   Matrix Market I/O example program
 *
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and copies it to stdout.  This porgram does nothing useful, but
 *   illustrates common usage of the Matrix Matrix I/O routines.
 *   (See http://math.nist.gov/MatrixMarket for details.)
 *
 *   Usage:  a.out [filename] > output
 *
 *
 *   NOTES:
 *
 *   1) Matrix Market files are always 1-based, i.e. the index of the first
 *      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
 *      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
 *      to files.
 *
 *   2) ANSI C requires one to use the "l" format modifier when reading
 *      double precision floating point numbers in scanf() and
 *      its variants.  For example, use "%lf", "%lg", or "%le"
 *      when reading doubles, otherwise errors will occur.

http://math.nist.gov/MatrixMarket/
http://math.nist.gov/MatrixMarket/mmio-c.html


gcc -O2 FIS_CG.c mmio.c operations.c -Wall -o csr.x -lm


./csr.x test.mtx
./csr.x tests.mtx
./csr.x orsirr_1.mtx
./csr.x s3rmt3m3.mtx
 */

/*
   gcc cg.c mmio.c operations.c -Wall -o cg.x -lm
   gcc -fopenmp -O2 cg.c mmio.c operations.c -Wall -o cg.x -lm


   ./csr.x test.mtx
   ./csr.x tests.mtx
   ./csr.x orsirr_1.mtx
   ./cg.x s3rmt3m3.mtx
   s3rmt3m3.mtx

	bcsstk01.mtx
	bcsstk12.mtx


   */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mmio.h"
#include "operations.h"
#include <time.h>
#include <omp.h>


double tol = 1e-08;



int main(int argc, char *argv[])
{
	int ret_code;
	MM_typecode matcode;
	FILE *f,*out, *out2;
	int M, N, nz, symm = 0;
	int i, j, m, *I, *J, *Ic, *Jc, *J_csr, *I_csc;
	double *val, *valc, *y, str_t, end_t, run_i;

str_t=omp_get_wtime();

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


	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
			mm_is_sparse(matcode) )
	{
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	/* find out size of sparse matrix .... */

	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
		exit(1);



	printf("Input check is finished. \n");
	/* reserve memory for matrices in CSR format */


	if (mm_is_general(matcode) == 1)
	{
		I = (int *) calloc(nz, sizeof(int));
		Ic = (int *) calloc(nz, sizeof(int));
		J_csr = (int *) calloc((M+1), sizeof(int));
		I_csc = (int *) calloc((M+1), sizeof(int));
		J = (int *) calloc(nz, sizeof(int));
		Jc = (int *) calloc(nz, sizeof(int));
		val = (double *) calloc(nz, sizeof(double));
		valc = (double *) calloc(nz, sizeof(double));
		y = (double *) calloc(N, sizeof(double));
	}

	if (mm_is_symmetric(matcode) == 1)
	{
		I = (int *) calloc(nz, sizeof(int));
		Ic = (int *) calloc(nz, sizeof(int));
		I_csc = (int *) calloc((M+1), sizeof(int));
		J_csr = (int *) calloc((M+1), sizeof(int));
		J = (int *) calloc(nz, sizeof(int));
		Jc = (int *) calloc(2*nz, sizeof(int));
		val = (double *) calloc(nz, sizeof(double));
		valc = (double *) calloc(2*nz, sizeof(double));
		y = (double *) calloc(N, sizeof(double));
		symm = 1;
	}

	// To make it sure, I set large size on below arrays.



	/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
	/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
	/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

	// I: X, J: Y reading data
	for (i=0; i<nz; i++)
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

	// if the matrix is symmetric, nz will be changed.
	printf("Converting COO into COO2. \n");
	nz = coo_to_csr(Ic, I, Jc, J_csr, J, valc, I_csc, val, M, nz, symm);
	printf("Finished converting COO2 into CSR. \n");


	/* Random numbers in a vector X*/


	/* Start GMRES method from now */


	/* The real GMRES codes starts from here*/

	// r_s and e has large size. That is because we don't know whether the number of iteration will be bigger or smaller than M.

	double *Y = (double *) calloc(M, sizeof(double)); // original equ's output vector
	double *r = (double *) calloc(M, sizeof(double)); // residual vector
	double *x = (double *) calloc(M, sizeof(double)); // initial guess & solution
	double *x_r = (double *) calloc(M, sizeof(double)); // Known real solution
	double *temp_arr = (double *) calloc(M, sizeof(double)); // temporary array
	double *temp_arr2 = (double *) calloc(M, sizeof(double)); // temporary array
	double *p = (double *) calloc(M, sizeof(double)); // search direction
	double alpha=0; // alpha values
	double norm_r2; // norm of r_m+1
	double r_s;
	double e;



	// Initializing initial guess, known solution and krylov space;
	for (i=0; i<N; i++)
	{
		x[i] = 0;
		x_r[i] = 1;
	}
	printf("Real solution and initial guess initialized \n");

	// temp_arr = A*x
		matrix_vector(J_csr, Jc, I_csc, I, val, valc, x, symm, temp_arr, M);

	// Initializing Y = Ax_r (know solution is x[i] = 1), Y is as good as b

		matrix_vector(J_csr, Jc, I_csc, I, val, valc, x_r, symm, Y, M);

		// residual is okay.
		printf("initial residual\n");
	for (i=0; i<M; i++)
	{
		r[i] = Y[i] - temp_arr[i];
//		printf("r[%d] = %e\n",i, r[i]);
	}

// measure runtime
run_i = clock();
	/* Iteration to do Conjugate Gradient Method starts now */
	double norm_r = norm(r, M);
	double r0 = norm_r;
	double beta;
	r_s = r0;
	printf("norm_r = %e \n",norm_r);
    for (i=0;i <M; i++)
    {
    	p[i] = r[i];
    }
    m=0;

	out=fopen("r_norm.txt","w");
//	out2=fopen("e_norm.txt","w");


    while ((r_s/r0 > tol)) // && ( m <=2)
    {
	// temp_arr2 = A*p
        matrix_vector(J_csr, Jc, I_csc, I, val, valc, p, symm, temp_arr2, M);
/*        for (i=0;i<M;i++)
        printf("temp_arr2[%d] = %e\n",i,temp_arr2[i]);
*/
        norm_r = dot_p(r,r,M); // Norm^2 is saved
        alpha = norm_r / dot_p(p,temp_arr2,M);
//        printf("dot_p(p,temp_arr2,M) %e\n",dot_p(p,temp_arr2,M) );
//        printf("alpha = %e\n",alpha);

//#pragma omp parallel for
        for (i=0;i <M; i++)
        {
        	x[i] = x[i] + alpha * p[i];
        }

        norm_r2 = 0;

       	for (i=0; i<M; i++)
        {
        	r[i] = r[i] - alpha*temp_arr2[i];
        	norm_r2 = norm_r2 + r[i]*r[i];
        }
//       	printf("norm_r2 = %e \n",norm_r2);
        beta= norm_r2/norm_r; // actually it is (rm+1,rm+1)/(rm,rm)

        printf("%d-th beta = %e\n",m,i,beta);

        for (i=0; i<M; i++)
        {
        	p[i] = r[i] + beta*p[i];
        }

        r_s = sqrt(norm_r2);
		fprintf(out,"%e \n", r_s);

//        e[m+1] = r_s[m+1]/r0;

        if (fabs(r_s/r0) < tol)
		{
        	m++;
        	printf("Solution converged!!\n");
        	break;
		}
        else
        {
        	printf("m = %d\n",m);
        m++;
        }
	}

    end_t=omp_get_wtime();
    fclose(out); // finish writing norm of r

    printf("steps to convergence %d \n",m);
    printf("runtime %e \n",end_t-str_t);


	/************************/
	/* now write out matrix in CSR format */
	/************************/
	printf("Finished running and writing files\n");

	out=fopen("J_csr.txt","w");
	fprintf(out, "Printing J_csr, N= %d\n\n", N);
	for (i=0; i < N+1; i++)
	{
		J_csr[i]++; // adjust from 0-based to 1-based.
		fprintf(out,"J_csr[%d] = %d \n", i, J_csr[i]);
	}
	fclose(out);

	printf("r0 = %e\n", norm_r);
	printf("Finished and getting into storing results\n");


	out=fopen("J_csr.txt","w");
	for (i=0; i < N+1; i++)
	{
		J_csr[i]++; // adjust from 0-based to 1-based.
		fprintf(out,"J_csr[%d] = %d \n", i, J_csr[i]);
	}
	fclose(out);




	out=fopen("J.txt","w");
	for (i=0; i < nz; i++)
	{
		Jc[i]++; // adjust from 0-based to 1-based.
		fprintf(out,"J[%d] = %d \n", i, Jc[i]);
	}
	fclose(out);


	out=fopen("val.txt","w");
	for (i=0; i < nz; i++)
	{
		fprintf(out," %e \n", valc[i]);
	}
	fclose(out);


	out=fopen("x.txt","w");
	for (i=0; i < M; i++)
	{
		fprintf(out,"%e \n", x[i]);
	}
	fclose(out);

	out=fopen("x_r.txt","w");
	for (i=0; i < M; i++)
	{
		fprintf(out,"%e \n", x_r[i]);
	}
	fclose(out);

	out=fopen("e.txt","w");
	for (i=0; i < M; i++)
	{
		fprintf(out,"%e \n", x_r[i]);
	}
	fclose(out);



	printf("Finished writing output files \n");

	free(I);
	free(J_csr);
	free(J);
	free(val);
	free(y);
	free(Ic);
	free(Jc);
	free(valc);


	free(Y);
	free(r);

// I got memory error on r_s
	printf("check point\n");

	free(x);
	free(x_r);
	free(temp_arr);
	free(temp_arr2);
	free(p);



printf("Congratulations, everything is well done\n");
	return 0;
}


