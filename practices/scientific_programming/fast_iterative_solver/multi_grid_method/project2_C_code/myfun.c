
double pi = 3.14159265359;
#include <stdio.h>
#include <string.h> // memcpy
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "myfun.h"
#include <unistd.h> // lseek,read,write
#include <fcntl.h> // open etc.
#include <err.h>

void Initialization(double* u_s, double* f, int N)
{
    double h = (1/(double)N);

    for (int i=1; i<N; i++)
    {
	for (int j=1; j<N; j++)
	{
	    double x = i*h;
	    double y = j*h;
	    u_s[i*(N+1) + j] = sin(2*pi*x)*sin(2*pi*y);
	    f[i*(N+1) + j] = 8*pi*pi*u_s[i*(N+1) + j];
	}
    }
}

void swapbytes(char* array, int nelem, int elemsize) {
    register int sizet, sizem, i, j;
    char* bytea;
    sizet = elemsize;
    sizem = sizet - 1;
    bytea = malloc(sizet);
    for(i = 0; i < nelem; ++i) {
	memcpy((void*) bytea, (void*) (array + i*sizet), sizet);
	for(j = 0; j < sizet; ++j) {
	    array[i*sizet + j] = bytea[sizem - j];
	}
    }
    free(bytea);
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


void M_print(double* A, int M, int N)
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

void scale(double* A, int M, double factor)
{
    int i;
    for (i=0; i<M;i++)
    {
	A[i] = A[i]*factor;
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
		M_print((B), M, M);

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
	M_print((result), M, M);

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
    printf("GS smoothing done iteration = %d norm = %e \n",iteration,  norm);

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


    GS(u, f, nu1, N);

    // Implement r_l = f- A_l*u_l (A_l: is Lagrangian operator)
    // u_t = A_l*u_t
    Laplacian(u_t, u, h, N);

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



