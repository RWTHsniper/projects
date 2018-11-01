

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "myfun.h"

int main(void)
{
    int l=7;
    int N= pow(2,l);
    int Nc = (N/2);
    double h = (1/(double)N);
    double (*u) = (double(*)) calloc((N+1)*(N+1), sizeof(double));
    // coarse mesh
    double (*u_c) = (double(*)) calloc((Nc+1)*(Nc+1), sizeof(double));
    // real solution
    double (*u_s) = (double(*)) calloc((N+1)*(N+1), sizeof(double));
    double (*f) = (double(*)) calloc((N+1)*(N+1), sizeof(double));

    printf("N = %d\n",N);
    printf("Nc = %d\n",Nc);
    printf("h = %e\n",h);

    // f and u_s initialization
    Initialization(u_s, f, N);

    int gamma = 2;
    int nu1 = 7, nu2 = 1, iteration = 30;

    double r0 = abs_max(f, (N+1)*(N+1));
    double r[1 + iteration];
    double r_time[1 + iteration];
    double sta_t=clock();

    r[0] = 1;

    // Perform multigrid method for a certain number of times.
    for (int m=0; m < iteration; m++)
    {
	MG(l, u, f, N, gamma, nu1, nu2);
	r[m+1] = r_inf_norm(u, f, N)/r0;
	r_time[m+1] = (clock()- sta_t)/CLOCKS_PER_SEC;
    }

    char aa[20] = "error";
    char bb[20] = "r_time";
    char cc[20] = ".txt";
    char nu1_c[10];
    sprintf(nu1_c,"%d",nu1);

    strcat(aa,nu1_c);
    strcat(aa,cc);

    strcat(bb,nu1_c);
    strcat(bb,cc);

    // save error and runtime as text file.
    M_fprint(aa, r, iteration+1, 1);
    M_fprint(bb, r_time, iteration+1, 1);

    if ((r[iteration+1]) < 1e-8)
    {
	printf("\nconverged! Relative error = %e \n\n", (r[iteration+1]/r0));
    }
    else
    {
	printf("\nFailed to converge, %e\n\n", (r[iteration+1]/r0));
    }


    free (u);
    free (u_c);
    free (u_s);
    free (f);

    return 0;
}



