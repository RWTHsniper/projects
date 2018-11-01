
double pi = 3.14159265359;
#include <stdio.h>
#include <string.h> // memcpy
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "myfun.h"
#include "mmio.h"
#include <unistd.h> // lseek,read,write
#include <fcntl.h> // open etc.
#include <err.h>


#include <iostream>
#include <sstream>
#include <string>

using namespace std;

void D_print(double* A, int M, int N)
{

    int i, j;
    for (i=0;i<M;i++)
    {
	for (j=0; j<N; j++)
	{
	    printf("%.3e ", A[i*N+j]);

	}
	printf("\n");
    }
}

void I_print(int* A, int M, int N)
{

    int i, j;
    for (i=0;i<M;i++)
    {
	for (j=0; j<N; j++)
	{
	    printf("%d ", A[i*N+j]);

	}
	printf("\n");
    }
}

void M_fprint(char* name, double* A, int M, int N)
{
    FILE* f;

    f=fopen(name,"w");

    int i, j;
    for (i=0;i<M;i++)
    {
	for (j=0; j<N; j++)
	{
	    fprintf(f, "%e ", A[i*N+j]);

	}
	fprintf(f,"\n");
    }
    fclose(f);

    /*
       Instruction how to put the name.	

       char aa[20] = "error";
       char bb[20] = "r_time";
       char cc[20] = ".txt";
       char nu1_c[10];
       sprintf(nu1_c,"%d",nu1);

       strcat(aa,nu1_c);
       strcat(aa,cc);

       strcat(bb,nu1_c);
       strcat(bb,cc);
       M_fprint(aa, r, iteration+1, 1);


       c++ version

       string aa = "error";
       string bb = "r_time";
       string cc = ".txt";
       string text = to_string(5);
       aa = aa+text+cc;
       M_fprint((char*)aa.c_str(), valc, M, 1);



     */


}


void M_I_fprint(char* name, int* A, int M, int N)
{
    FILE* f;

    f=fopen(name,"w");

    int i, j;
    for (i=0;i<M;i++)
    {
	for (j=0; j<N; j++)
	{
	    fprintf(f, "%d ", A[i*N+j]);

	}
	fprintf(f,"\n");
    }
    fclose(f);

    /*
       Instruction how to put the name.	

       char aa[20] = "error";
       char bb[20] = "r_time";
       char cc[20] = ".txt";
       char nu1_c[10];
       sprintf(nu1_c,"%d",nu1);

       strcat(aa,nu1_c);
       strcat(aa,cc);

       strcat(bb,nu1_c);
       strcat(bb,cc);
       M_fprint(aa, r, iteration+1, 1);


       c++ version

       string aa = "error";
       string bb = "r_time";
       string cc = ".txt";
       string text = to_string(5);
       aa = aa+text+cc;
       M_fprint((char*)aa.c_str(), valc, M, 1);



     */


}


double abs_max(double *A, int M) // find abs_max
{
    double max=0;
    for (int i=0; i<M; i++)
    {
	if (max < fabs(A[i]))
	{
	    max = fabs(A[i]);
	}

    }
    return max;
}

double norm(double* A, int M)
{
    double sum = 0;
    for (int i=0; i<M; i++)
    {
	sum = sum+ A[i]*A[i];
    }
    sum= pow(sum,0.5);

    return sum;
}

void Unit_v(double* A, double* B, int M)
{//	output unit_v, input_v, size of vector
    double n = norm(B,M);

    for (int i=0; i<M;i++)
    {
	A[i] = B[i]/n;
    }

}


double M_inf_norm(double* A, int r, int c)
{
    double sum, temp;
    int i, j;

    sum = 0;
    for (i=0; i<r; i++)
    {
	temp = 0;
	for (j=0;j<c; j++)
	{
	    temp = temp + fabs(A[c*i+ j]);
	}

	if (temp > sum )
	{
	    sum = temp;
	}
    }

    return sum;
}

void scale(double* A, double* B, int M, double factor)
{
    int i;
    for (i=0; i<M;i++)
    {
	A[i] = B[i]*factor;
    }
}


void copy(double* B, double*A, int M)
{
    int i;
    for (i=0; i<M;i++)
    {
	B[i] = A[i];
    }
}

void exchange(double *A, double *B, int M )
{
    double (*C) = (double(*)) calloc(M, sizeof(double));

    for (int i=0; i<M; i++)
    {
	C[i] = A[i];
	A[i] = B[i];
	B[i] = C[i];
    }


    free(C);
}

double dot_p(double *A, double *B, int M)
{
    double res=0;
    for (int i=0;i<M;i++)
    {
	res = res + A[i]*B[i];
    }
    return res;
}



void M_transpose(double *B, double *A, int M)
{
    double (*C) = (double(*)) calloc(M*M, sizeof(double));


    for (int i=0; i<M; i++)
	for (int j=0; j<M; j++)
	{
	    C[M*j+i] = A[M*i+j];
	}
    copy(B, C, M*M);

    free (C);
}


void add(double* C, double* A, double* B, int M)
{
    // A + B = C
    int i;
    for (i=0; i<M;i++)
    {
	C[i] = A[i]+B[i];
    }
}

void subtract(double* C, double* A, double* B, int M)
{
    // A - B = C
    int i;
    for (i=0; i<M;i++)
    {
	C[i] = A[i]-B[i];
    }
}




void Inv_sign(double*A, int M)
{
    //  A is array, M is total size of terms to be *-1
    for (int i=0; i<M; i++)
    {
	A[i] = -A[i];
    }
}

void swapbytes(char* array, int nelem, int elemsize) {
    register int sizet, sizem, i, j;
    char* bytea;
    sizet = elemsize;
    sizem = sizet - 1;
    //    bytea = malloc(sizet);
    bytea = (char*)malloc(sizet);
    for(i = 0; i < nelem; ++i) {
	memcpy((void*) bytea, (void*) (array + i*sizet), sizet);
	for(j = 0; j < sizet; ++j) {
	    array[i*sizet + j] = bytea[sizem - j];
	}
    }
    free(bytea);
}

void address_swaper(int **A, int **B)
{
    int *temp = *A;
    *A = *B;
    *B= temp;
}

void address_swaper(double **A, double **B)
{
    double *temp = *A;
    *A = *B;
    *B= temp;
}


// Use this routine in order to read a file containing double data
void read_data(
	double* data,     // Output: Array containing the data
	size_t nn,        // Input: Number of nodes/elements to read
	unsigned entries, // Input: Number of data entries (DOFs or coordinates) per node/element
	char* file_name   // Input: Path to the data file
	){
    int fd;
    fd = open(file_name, O_RDONLY);
    if (fd == -1) {
	err(1,"Opening %s for reading", file_name);
    }

    size_t bytestoread = nn*entries*sizeof(double);// CODE HERE: Number of bytes to read
    size_t bytesread = read(fd, data, bytestoread);
    if (bytesread == -1) {
	err(1, "%s: Reading data failed", file_name);
    } else if (bytesread != bytestoread) {
	printf("%s: Only %zd bytes of %zu bytes read!\n", file_name, bytesread, bytestoread);
	exit(1);
    }

    close(fd);

    swapbytes((char*)data, nn*entries, sizeof(double));

}

// Use this routine in order to write the projected data
void write_data(
	double* data,     // Input: Data to be written
	size_t nn,        // Input: Number of nodes to write
	unsigned entries, // Input: Number of degrees of freedom
	char* file_name   // Input: Path to the data file
	){

    swapbytes((char*)data, nn*entries, sizeof(double));
    int fd;
    fd = open(file_name, O_RDWR | O_CREAT, 0666);
    if (fd == -1) {
	err(1,"Opening %s for writing", file_name);
    }

    size_t bytestowrite = nn*entries*sizeof(double);// CODE HERE: Number of bytes to read
    size_t byteswritten = write(fd, data, bytestowrite);
    if (byteswritten == -1) {
	err(1, "%s: Writing data failed", file_name);
    } else if (byteswritten != bytestowrite) {
	printf("%s: Only %zd bytes of %zu bytes written!\n", file_name, byteswritten, bytestowrite);
	exit(1);
    }

    close(fd);
}


void Init(double*A, int M, double N) 	// Initialize arrays
{/* A: array, M: total number of elements to be initialized
N: initialized number */
    for(int i=0; i<M; i++)
    {
	A[i] = N;
    }
}



double determinant(double* A, int M)
{
    int n = M;
    double tol = 1e-20;
    double det=1;
    double (*B) = (double(*)) calloc(M*M, sizeof(double));
    if(n==1) return A[0*M + 0];
    if(n==2) return (A[0*M + 0]*A[1*M + 1] - A[0*M + 1]*A[1*M + 0]);
    if(n==3) return (A[0*M + 0]*((A[1*M + 1]*A[2*M+2]) - (A[2*M+1]*A[1*M+2])) -A[0*M+1]*(A[1*M+0]*A[2*M+2] - A[2*M+0]*A[1*M+2]) + A[0*M+2]*(A[1*M+0]*A[2*M+1] - A[2*M+0]*A[1*M+1]));


    // Do row reduction to make echelon form
    copy(B, A, M*M); // Copy mat onto B


    //printf("trying to make echelon form\n");
    for (int k=0; k<M-1; k++) // k: reference row
    {
	for (int i=k+1; i<M; i++) // i: target row
	{

	    if (fabs(B[k*M+k]) < tol)
	    {
		//printf("[%d %d]exchaning rows\n",i,k);
		for (int l=k+1; l<M; l++)
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

		printf("singular pivot, i, k, M, B1, B2 = %d %d %d %e, %e\n",i, k, M, B[i*M+k],B[k*M+k]);

		printf("B printing \n\n");
		D_print((B), M, M);

		return 0;
	    }

	    double fact = B[i*M+k]/B[k*M+k];
	    //printf("[%d, %d] %e\n",i, k, fact);
	    for (int j=k; j<M; j++) // calculating jth column
	    {
		B[i*M+j] = B[i*M+j] -fact*B[k*M+j];
	    }


	}
    }

    for (int i=0; i<M; i++)
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

void check_value(double *A, int M) // checks every terms in an array
{

    for (int i=0; i<M; i++)
    {
	if (isfinite(A[i]) == 0)
	{
	    printf("Check value fucked up.");
	    printf("A[%d] = %e",i,A[i]);
	    exit(-1);
	}
    }
}


void D_matrix_vector(double* y, double* A, int M, int N, double* x)
{ //                 	output inp Matrix    row   column   inp      
    int i, j;
    double tol = 1e-8;

    // Initializing y where the result will be installed;
    for (i=0; i<M; i++)
    {
	y[i] = 0;
    }

    // running program
    for (i=0;i<M; i++)
    {
	for (j=0; j<N; j++)
	{
	    if (fabs(x[j]) > tol)
	    {
		y[i] = y[i] + A[N*i+ j]*x[j];
	    }
	}
    }

    printf("\nM = %d \n",M);
    // output
    /*
       printf("\ny\n");
       for (j=0; j<M; j++)
       printf("%e\n", y[j]);
     */
}


double Inv_M(double *result, double *A, int M) // result storage, input matrix, size of matrix
    // input and output should be different
{
    printf("Inversion procedure start!\n");

    // M is dimension of input matrix A

    double *C = (double(*)) calloc((M-1)*(M-1), sizeof(double));
    double tol = 1e-16;
    double det=0;


    det = determinant(A, M);

    if (fabs(det) < tol)
    {
	printf("\nsingular matrix. Fucked up.\n");
	printf("tolerance = %e\n",tol);
	exit(-1);
    }

    // Build matrix of minors
    for (int i=0; i<M; i++)
    {

	for (int j=0; j<M; j++)
	{
	    int count=0;
	    // cofactor amtrix for B[i][j]
	    for (int k=0; k<M; k++)
	    {
		for (int l=0; l<M; l++)
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

    M_transpose(result, result, M);


    free (C);


    double inv_check= determinant(result,M);
    if ((1 - tol >det*inv_check) || (1 + tol < det*inv_check))
    {
	printf("\n[%e %e]Inverse matrix calculation error\n", det, inv_check);
	printf("tolerance = %e\n",tol);
	printf("result printing \n\n");
	D_print((result), M, M);

    }

    return det;
}





void GS(double* u0, double* f, int nu, int N)
{
    int iteration=0, i, j, k, l, m;
    double h = (1/(double)N);
    double tol = 1e-10;
    double norm = 1;
    double (*u_t) = (double(*)) calloc((N+1)*(N+1), sizeof(double));

    /*
       printf("u0 printing \n\n");
       D_print((u0), N+1, N+1);
     */

    for (m=0; m<nu; m++)
    {iteration++;
	// Load values onto u_t
	copy(u_t, u0, (N+1)*(N+1));

	for (j=1; j<=N-1;j++)
	{
	    for (i=1; i<=N-1;i++)
	    {
		u0[(N+1)*i + j] = (h*h*f[(N+1)*i + j] + u0[(N+1)*(i-1) + j] + u0[(N+1)*i + j-1] + u0[(N+1)*(i+1) + j] + u0[(N+1)*i + j+1])/4;

	    }
	}
    }
    //	printf("GS smoothing done iteration = %d norm = %e \n",iteration,  norm);

    free(u_t);
}




void Restriction(double* u_c, double* u, int N)
{
    int ii, jj, Nc = N/2;

    for (int i=1; i<=Nc-1;i++)
    {
	ii = 2*i;
	for (int j=1; j<=Nc-1;j++)
	{
	    jj = 2*j;
	    u_c[(Nc+1)*i + j] = (u[(N+1)*(ii-1) + (jj-1)] +2*u[(N+1)*(ii) + (jj-1)] + u[(N+1)*(ii+1) + (jj-1)] + 2*u[(N+1)*(ii-1) + (jj)] +4*u[(N+1)*(ii) + (jj)] + 2*u[(N+1)*(ii+1) + (jj)]
		    + u[(N+1)*(ii-1) + (jj+1)] + 2*u[(N+1)*(ii) + (jj+1)] + u[(N+1)*(ii+1) + (jj+1)])/16;
	}
    }
}



void Prolongation(double* u, double* u_c, int N)
{
    int Nc = N/2;
    // initialize u
    for (int i=1; i<=N-1; i++)
    {
	for (int j=1; j<=N-1; j++)
	{
	    u[(N+1)*i + j] = 0;
	}
    }

    for (int i=1; i<=Nc-1; i++)
    {
	int ii = 2*i;
	for (int j=1; j<=Nc-1; j++)
	{
	    int jj = 2*j;
	    double temp = (u_c[(Nc+1)*(i)+ (j)]);
	    u[(N+1)*(ii-1)+ (jj-1)] = u[(N+1)*(ii-1)+ (jj-1)] + (double)1/4*temp;
	    u[(N+1)*(ii-1)+ (jj+1)] = u[(N+1)*(ii-1)+ (jj+1)] + (double)1/4*temp;
	    u[(N+1)*(ii+1)+ (jj+1)] = u[(N+1)*(ii+1)+ (jj+1)] + (double)1/4*temp;
	    u[(N+1)*(ii+1)+ (jj-1)] = u[(N+1)*(ii+1)+ (jj-1)] + (double)1/4*temp;

	    u[(N+1)*(ii)+ (jj-1)] = u[(N+1)*(ii)+ (jj-1)] + (double)1/2*temp;
	    u[(N+1)*(ii)+ (jj+1)] = u[(N+1)*(ii)+ (jj+1)] + (double)1/2*temp;
	    u[(N+1)*(ii-1)+ (jj)] = u[(N+1)*(ii-1)+ (jj)] + (double)1/2*temp;
	    u[(N+1)*(ii+1)+ (jj)] = u[(N+1)*(ii+1)+ (jj)] + (double)1/2*temp;

	    u[(N+1)*(ii)+ (jj)] += temp;

	}
    }
}

void Laplacian(double* u_t, double* u, double h, int N)
{
    //Init(u_t, (N+1)*(N+1), 0); 	// Initialize arrays, size of array, initializing value

    //	for (int i=0; i<N+1; i++)
    for (int i=1; i<N; i++)
    {
	//		for (int j=0; j<N+1; j++)
	for (int j=1; j<N; j++)
	{
	    u_t[(N+1)*(i) + (j)] =	-(u[((N+1)*(i-1) + j)] -4 * u[((N+1)*(i) + j)] + u[((N+1)*(i+1) + j)] + u[((N+1)*(i) + (j-1))] + u[((N+1)*(i) + (j+1))]);
	    u_t[(N+1)*(i) + (j)] = u_t[(N+1)*(i) + (j)]/(h*h);

	}
    }

}


double r_inf_norm(double *u, double* f, int N)
{
    double h = (1/(double)N);
    double max = 0;
    double (*Lu) = (double(*)) calloc((N+1)*(N+1), sizeof(double));

    Laplacian(Lu, u, h, N);
    subtract(Lu, f, Lu, (N+1)*(N+1)); // A -B = C, M:length of array

    max = abs_max(Lu, (N+1)*(N+1)); // find abs_max

    free (Lu);
    return max;
}


void Inv_Laplacian(double *result, double *u, double h, int N)
{
    printf("Inv laplacian start\n");


    // build up laplacian operator
    int size  = (N+1)*(N+1);
    double (*L) = (double(*)) calloc(size*size, sizeof(double));
    double (*L_inv) = (double(*)) calloc(size*size, sizeof(double));

    L[size*0+0] = (double)4/(h*h);
    L[size*0+1] = (double)-1/(h*h);
    L[size*0+(N+1)] = (double)-1/(h*h);

    for (int i=1; i<size-1; i++)
    {
	L[size*i+i] = (double)4/(h*h);
	L[size*i + i+1] = (double)-1/(h*h);
	L[size*i + i-1] = (double)-1/(h*h);
	if ((i+(N+1)) <= size-1)
	{
	    L[size*i + i+(N+1)] = (double)-1/(h*h);
	}

	if ((i-(N+1)) >= 0)
	{
	    L[size*i + i-(N+1)] = (double)-1/(h*h);
	}
    }
    L[size*(size-1)+size-1] = (double)4/(h*h);
    L[size*(size-1) + size-2] = (double)-1/(h*h);
    L[size*(size-1) + size-1 -(N+1)] = (double)-1/(h*h);

    Inv_M( L_inv, L, size);

    double (*u_t) = (double(*)) calloc(size, sizeof(double));
    copy(u_t, u, size); // Copy u onto u_t

    D_matrix_vector(result, L_inv, size, size, u_t);

    free (L);
    free (L_inv);
    free (u_t);
}


// I think f, u, u_c should be depend on size of mesh
void MG(int l, double *u, double *f, int N, int gamma, int nu1, int nu2)
{// l: level, u: grid, f, N: N at l-th level, gamma: #iteration at each level
    int Nc = N/2;
    double *u_c = (double(*)) calloc((Nc+1)*(Nc+1), sizeof(double));
    double (*r_l) = (double(*)) calloc((N+1)*(N+1), sizeof(double));
    double (*r_l2) = (double(*)) calloc((Nc+1)*(Nc+1), sizeof(double));
    double (*e_l) = (double(*)) calloc((N+1)*(N+1), sizeof(double));
    double (*e_l2) = (double(*)) calloc((Nc+1)*(Nc+1), sizeof(double));
    double (*u_t) = (double(*)) calloc((N+1)*(N+1), sizeof(double));
    double h = (1/(double)N);
    double h2 = (double)h/2;

    printf("\nstart %d-th level, h = %e N = %d \n",l,h,N);

    /*
       printf("f printing \n\n");
       M_print((f), N+1, N+1);
     */

    GS(u, f, nu1, N);

    /*	printf("u printing \n\n");
	M_print((u), N+1, N+1);
     */
    /*
       printf("\nu_t printing before Laplacian \n\n");
       M_print((u_t), N+1, N+1);
     */

    // Implement r_l = f- A_l*u_l (A_l: is Lagrangian operator)
    // u_t = A_l*u_t
    Laplacian(u_t, u, h, N);
    /*
       printf("u_t printing after Laplacian \n\n");
       M_print((u_t), N+1, N+1);
     */

    subtract(r_l, f, u_t, (N+1)*(N+1));
    // r_l2 = Restriction r_l
    Restriction(r_l2, r_l, N);

    if (l==1)
    {
	printf("reached at just before the lowest level Nc = %d\n", Nc);
	GS(e_l2, r_l2, 1, Nc);
	//Inv_Laplacian( e_l2, r_l2, h2, Nc);
	// Jinxuan told me that Gauss relaxation is inverse opeation of laplacian
	// e_l2 = -e_l2
	Inv_sign( e_l2, (Nc+1)*(Nc+1));
	printf("The lowest error equation is solved!\n");
    }
    else
    {
	Inv_sign( r_l2, (Nc+1)*(Nc+1));
	Init(e_l2, (Nc+1)*(Nc+1), 0);
	for (int j=0; j< gamma; j++)
	{
	    MG(l-1, e_l2, r_l2, Nc, gamma, nu1, nu2);
	}
    }
    // Continue writing from here
    Prolongation(e_l, e_l2, N);

    subtract(u, u, e_l, (N+1)*(N+1));

    GS(u, f, nu2, N);


    free (u_c);
    free (r_l);
    free (r_l2);
    free (e_l);
    free (e_l2);
    free (u_t);
}



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
    // I_csc, J_csr will be built
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


void CSR_matrix_vector(int* J_csr, int* Jc, int* I_csc, int* I, double* val, double* valc, double* Xv, int symm, double* y, int M)
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

// Unfinished
void CSR_vector_matrix(int* J_csr, int* Jc, int* I_csc, int* I, double* val, double* valc, double* Xv, int symm, double* y, int M)
{
    //                CSR form, Jc: column, value vector, Xv: input vector, y: output, M: number of rows
    int i, j;

    // Initializing y where the result will be installed;
    for (i=0; i<M; i++)
    {
	y[i] = 0;
    }

    if (symm ==0)
    {	// i-th column, j-th row
	for (i=0; i<M; i++)
	    for (j=I_csc[i]; j<I_csc[i+1]; j++)
	    {

		y[i] = y[i] + valc[ j ]*Xv[ Jc[j] ];
	    }
    }
    else
    {
	// do calculation in column using csc
	for (i=0; i<M; i++)
	    for (j=I_csc[i]; j<I_csc[i+1]; j++)
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
	// do calculation in row using csr

	for (i=0; i<M; i++)
	{
	    for (j=J_csr[i] + 1; j<J_csr[i+1]; j++)
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

// Unfinished
void CSR_print(int *J_csr, int *Jc, double* valc, int M, int N)
{

    for (int i=0; i<M; i++)
    {
	for (int j=J_csr[i] + 1; j<J_csr[i+1]; j++)
	{
	}
    }
}




void CSR_CG_solver(double *x, int *J_csr, int *Jc, int *I_csc, int *I, double *val, double *valc, int symm, double *Y, int M, double tol)
{		// x: solution, .........; Y: rhs of initial equation
    int i;

    double *temp_arr = (double *) calloc(M, sizeof(double)); // temporary array
    double *temp_arr2 = (double *) calloc(M, sizeof(double)); // temporary array
    double *p = (double *) calloc(M, sizeof(double)); // search direction
    double *r = (double *) calloc(M, sizeof(double)); // residual vector
    double alpha=0; // alpha values
    double norm_r2; // norm of r_m+1
    double r_s;
    double e;
    int m=0;

    // Initializing initial guess, known solution and krylov space;
    for (i=0; i<M; i++)
    {
	x[i] = 0;
    }

    // temp_arr = A*x
    CSR_matrix_vector(J_csr, Jc, I_csc, I, val, valc, x, symm, temp_arr, M);


    // printf("initial residual\n");
    subtract(r, Y, temp_arr, M); 

    double norm_r = norm(r, M);
    double r0 = norm_r;
    double beta;
    r_s = r0;

    copy(p, r, M); // Copy r onto p




    while ((r_s/r0 > tol)) // && ( m <=2)
    {
	// temp_arr2 = A*p
	CSR_matrix_vector(J_csr, Jc, I_csc, I, val, valc, p, symm, temp_arr2, M);
	/*        for (i=0;i<M;i++)
		  printf("temp_arr2[%d] = %e\n",i,temp_arr2[i]);
	 */
	norm_r = dot_p(r,r,M); // Norm^2 is saved
	alpha = norm_r / dot_p(p,temp_arr2,M);
	//        printf("dot_p(p,temp_arr2,M) %e\n",dot_p(p,temp_arr2,M) );
	//        printf("alpha = %e\n",alpha);

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

	for (i=0; i<M; i++)
	{
	    p[i] = r[i] + beta*p[i];
	}

	r_s = sqrt(norm_r2);

	if (fabs(r_s/r0) < tol)
	{
	    m++;
	    //	    printf("Solution converged!!\n");
	    break;
	}
	else
	{
	    m++;
	}
    }

    free(temp_arr);
    free(temp_arr2);
    free(p);
    free(r);
}

double power_iteration(int *I_csc, int *I, int *J_csr, int *Jc, double *val, double *valc, int symm, int M, int N, double tol)
{

    int k=1;
    FILE *f;
    double *z = (double *) calloc(M, sizeof(double));
    double *q = (double *) calloc(M, sizeof(double));
    double *temp = (double *) calloc(M, sizeof(double));
    double error = 1;
    double time;
    double lam=0, prev_lam=0;
    int write=1;

	// initialize q
    for (int i=0; i<M;i++)
    {
	q[i] = 1/sqrt((double)M);
    }


    if (write ==1)
    {
	string a = "power";
	string b = to_string(M);
	a = a + b;
	f=fopen((char*)a.c_str(),"w");
	time = clock();
    }
    while (error > tol)
	//	while (k <= 100)
    {
	// z_k=A*q_(k-1)	
	CSR_matrix_vector(J_csr, Jc, I_csc, I, val, valc, q, symm, z, M);

	// q = z/||z||
	Unit_v(q, z, M);

	// temp = Aq
	CSR_matrix_vector(J_csr, Jc, I_csc, I, val, valc, q, symm, temp, M);

	// lam = dot(q,(Aq)) = q*(Aq)
	lam = dot_p(temp, q, M);
	error=fabs(lam-prev_lam);
	// runtime & error check.
	if (write ==1)
	{
	    fprintf(f,"%e %e\n",(clock()-time)/CLOCKS_PER_SEC,error);
	}
	if ( error< tol)
	{
	    printf("converged at %d\n",k);
	    break;
	}
	else
	{
	    prev_lam = lam;
	    k++;
	}
    }
    if (write ==1)
    {
	fclose(f);
    }
    printf("lambda = %e\n",lam);

    free(z);	
    free(q);
    free(temp);

    return lam;
}


double inv_power_iteration(int *I_csc, int *I, int *J_csr, int *Jc, double *val, double *valc, int symm, int M, int N, double tol)
{
    printf("\nLet's go for inverse method!\n");

    int k = 1;
    FILE *f;
	double time;
    double *z = (double *) calloc(M, sizeof(double));
    double *q = (double *) calloc(M, sizeof(double));
    double *temp = (double *) calloc(M, sizeof(double));
    double error = 1;
    double lam=0, prev_lam=0;
    int write=0;

    // initialize q
    for (int i=0; i<M;i++)
    {
	q[i] = 1/sqrt((double)M);
    }

    if (write == 0)
    {
	string a = "inv_power";
	string b = to_string(M);
	a = a + b;
	f=fopen((char*)a.c_str(),"w");
	time = clock();
    }

    while (error >= tol)
    {
	// z_k=(A^-1)*q_(k-1), A*z = q	
	CSR_CG_solver(z, J_csr, Jc, I_csc, I, val, valc, symm, q, M, 1e-8);


	Unit_v(q, z, M);

	// lambda = q'*(A^-1q)
	// A^-1q: A*temp = q
	CSR_CG_solver(temp, J_csr, Jc, I_csc, I, val, valc, symm, q, M, 1e-8);
	// q*(Aq)
	lam = dot_p(temp, q, M);
	error=fabs(lam-prev_lam);
	if (write == 0)
	{
	    fprintf(f,"%e %e\n",(clock()-time)/CLOCKS_PER_SEC,error);
	}
	if ( error< tol)
	{
	    printf("inverse power converged at %d\n",k);
	    break;
	}
	else
	{
	    prev_lam = lam;
	    k++;
	}
    }
    if (write == 0)
    {
	fclose(f);
    }
    lam = 1/lam;
    printf("lambda = %e\n",lam);


    free(z);	
    free(q);
    free(temp);

    return lam;

}

double CSR_lanczos(int *I_csc, int *I, int *J_csr, int *Jc, double *val, double *valc, int symm, int M, int N)
{
    double *v_p = (double *) calloc(M, sizeof(double));	// previour vector
    double *v = (double *) calloc(M, sizeof(double));	// current vector
    double *w = (double *) calloc(M, sizeof(double));
    double *beta = (double *) calloc(M+1, sizeof(double));
    double *alpha = (double *) calloc(M+1, sizeof(double));
    int write =1;

    // initialize v0
    for (int i=0; i<M;i++)
    {
	v[i] = 1/sqrt((double)M);
    }
    //	size of Krylov space for Lanczos method m: 30 50 75 100
    int m= 75;

    // I think I have problem with scale, add, dot_p operations.
    for (int j=0; j<m; j++)
    {
	// w = A * v
	CSR_matrix_vector(J_csr, Jc, I_csc, I, val, valc, v, symm, w, M);

	// I just use v_p = beta_j-1 * v_j-1 to save memory
	if (j > 0)
	{
	    // when j=0, beta = 0 and those operations are unnecessary.
	    scale(v_p, v_p, M, -beta[j-1]);
	    add(w, w, v_p, M);
	}
	alpha[j] = dot_p(v, w, M);
	// w = w - alph_j*vj
	// v_p = v_j-1 = -alpha_j*vj to save memory
	scale(v_p, v, M, -alpha[j]);
	add(w, w, v_p, M);
	beta[j]= norm(w, M);

	// let v_p to point v's data by simply changing address of pointer.
	address_swaper(&(v_p), &(v)); // time-saving than copying entire data.
	// update v = new v
	scale(v, w, M, (double)1/beta[j]);
    }

    if (write ==1)
    {
	M_fprint("alpha.txt", &(alpha[0]), m, 1);
	M_fprint("beta.txt", &(beta[0]), m-1, 1);
    }

    int T_nz = 1 + 2*(m-1);
    int T_symm = 1;
    int *T_I = (int *) calloc(T_nz, sizeof(int));
    int *T_I_csc = (int *) calloc(M+1, sizeof(int));
    int *T_J = (int *) calloc(T_nz, sizeof(int));
    int *T_J_csr = (int *) calloc(M+1, sizeof(int));
    double *T_I_val = (double *) calloc(T_nz, sizeof(double));
    double *T_J_val = (double *) calloc(T_nz, sizeof(double));


    // Express tridiagonal matrix using CSC/CSR

    // CSC converter

    T_I[0] = 0;
    T_I[1] = 1;
    T_I[T_nz-1] = m-1;

    T_I_val[0]=alpha[0];
    T_I_val[1]=beta[0];
    T_I_val[T_nz-1]=alpha[m-1];
    int count=2;
    for (int i=1; i<m-1; i++)
    {
	// diagonal term
	T_I[count]=i;
	T_I_val[count]= alpha[i];	
	count++;
	T_I[count]=i+1;
	T_I_val[count]= beta[i];	
	count++;
	T_I_csc[i]=2*i;
    }

    T_I_csc[0] = 0;
    T_I_csc[m-1] = T_nz-1;
    T_I_csc[m] = T_nz;


    // CSR converter

    T_J[0] = 0;
    T_J[T_nz-1] = m-1;

    T_J_val[0]=alpha[0];
    T_J_val[T_nz-1]=beta[m-2];
    count=1;
    for (int i=1; i<m; i++)
    {
	// off-diagonal term
	T_J[count]=i-1;
	T_J_val[count]= beta[i-1];	
	count++;
	// diagonal term
	T_J[count]=i;
	T_J_val[count]= alpha[i];	
	count++;
	T_J_csr[i]=2*i-1;
    }

    T_J_csr[0] = 0;
    T_J_csr[m] = T_nz;

    printf("Power iteration runs after Lanczos method\n");
    double tol=0;
    // tolerance modification according to Krylov space size.
    // next time, I am going to change it. Apply linear interpolation in each region.
    // y = y1 + (x-x1)/(x2-x1)*(y2-y1)
    if ((30<=m) && (m<50))
    {
	tol = 1e-2 + (m-30)/(50-30)*(1e-4 - 1e-2);
    }
    else if ((50<=m) && (m<75))
    {
	tol = 1e-4 + (m-50)/(75-50)*(1e-6 - 1e-4);
    }
    else if ((75<=m) && (m<100))
    {
	tol = 1e-6 + (m-75)/(100-75)*(1e-10 - 1e-6);	
    }
    else if (m>=100)
    {
	tol = 1e-10;
    }
    printf("[Lanczos] m = %d tol = %e\n",m,tol);

    double ret = power_iteration(T_I_csc, T_I, T_J_csr, T_J, T_I_val, T_J_val, symm, m, m, tol);

    free(v_p);
    free(v);
    free(w);
    free(beta);
    free(alpha);

    free(T_I);
    free(T_I_csc);
    free(T_J);
    free(T_J_csr);
    free(T_I_val);
    free(T_J_val);

    return ret;

}
