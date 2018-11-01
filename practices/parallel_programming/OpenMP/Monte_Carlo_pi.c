#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define TN 8
#define DEBUG 0
#define THRESHOLD (5e-3)

void print_matrix( double *phi, int n )
{
    int i, j;
    for ( i = 0; i < n; i++ )
    {
        printf("[" );
        for ( j = 0; j < n; j++ )
            printf(" %10.6f ", phi[j*n + i] );
        printf("]\n" );
    }
}

void init_phi( double * phi, int n )
{
    int i, j;

    // Interior points initialized to 50 degrees
    for ( i = 1; i < n-1; i++ )
        for ( j = 1; j < n-1; j++ )
            phi[j*n+i] = 50.0;

    // Top, left, and right boundaries fixed at 100 degrees
    for ( i = 0; i < n; i++ )
    {
        phi[    0*n +    i] = 100.0;
        phi[(n-1)*n +    i] = 100.0;
        phi[    i*n +    0] = 100.0;
    }
    // Bottom boundary fixed at 0 degrees
    for ( i = 0; i < n; i++ )
        phi[    i*n +(n-1)] = 0.0;
}


double compute_seq( int n, int *niters )
{
    double *phi_cur, *phi_next, *tmp;
    double conv, phi;
    int i, j;

    phi_cur  = (double *) malloc ( n * n * sizeof(double) );
    phi_next = (double *) malloc ( n * n * sizeof(double) );

    init_phi(phi_cur, n);
    init_phi(phi_next, n);

    *niters = 0;
    while( 1 )
    {
        (*niters)++;
#if DEBUG
        printf("Iteration %d\n", *niters);
        print_matrix( phi_cur, n );
        sleep(1);
#endif
        // Compute next (new) phi from current (old) phi
omp_set_num_threads(TN);

{

        for( j = 1; j < n-1; j++ )
            for( i = 1; i < n-1; i++ )
                phi_next[ j*n+i ] = ( phi_cur[ (j-1)*n+i ] +
                                      phi_cur[ j*n+(i-1) ] +
                                      phi_cur[ j*n+(i+1) ] +
                                      phi_cur[ (j+1)*n+i ] ) / 4;
}
        // If converged, we are done
        conv = 1;
        for( j = 1; j < n-1; j++ )
            for( i = 1; i < n-1; i++ )
                if( fabs( phi_next[ j*n + i ] - phi_cur[ j*n + i ] ) > THRESHOLD )
                    conv = 0;
        if ( conv )
            break;
        // Otherwise, swap pointers and continue
        tmp = phi_cur;
        phi_cur = phi_next;
        phi_next = tmp;
    }

    free( phi_cur );
    phi=*phi_next;
    free( phi_next );
return phi;
}


double compute_outer( int n, int *niters )
{
    double *phi_cur, *phi_next, *tmp;
    double conv, phi;
    int i, j;

    phi_cur  = (double *) malloc ( n * n * sizeof(double) );
    phi_next = (double *) malloc ( n * n * sizeof(double) );

    init_phi(phi_cur, n);
    init_phi(phi_next, n);

    *niters = 0;
    while( 1 )
    {
        (*niters)++;
#if DEBUG
        printf("Iteration %d\n", *niters);
        print_matrix( phi_cur, n );
        sleep(1);
#endif
        // Compute next (new) phi from current (old) phi
omp_set_num_threads(TN);
#pragma omp parallel
{
#pragma omp for private(i)
        for( j = 1; j < n-1; j++ )
            for( i = 1; i < n-1; i++ )
                phi_next[ j*n+i ] = ( phi_cur[ (j-1)*n+i ] +
                                      phi_cur[ j*n+(i-1) ] +
                                      phi_cur[ j*n+(i+1) ] +
                                      phi_cur[ (j+1)*n+i ] ) / 4;
}
        // If converged, we are done
        conv = 1;
        for( j = 1; j < n-1; j++ )
            for( i = 1; i < n-1; i++ )
                if( fabs( phi_next[ j*n + i ] - phi_cur[ j*n + i ] ) > THRESHOLD )
                    conv = 0;
        if ( conv )
            break;
        // Otherwise, swap pointers and continue
        tmp = phi_cur;
        phi_cur = phi_next;
        phi_next = tmp;
    }

    free( phi_cur );
    phi=*phi_next;
    free( phi_next );
return phi;
}


void compute_inner( int n, int *niters )
{
    double *phi_cur, *phi_next, *tmp;
    double conv;
    int i, j;

    phi_cur  = (double *) malloc ( n * n * sizeof(double) );
    phi_next = (double *) malloc ( n * n * sizeof(double) );

    init_phi(phi_cur, n);
    init_phi(phi_next, n);

    *niters = 0;
    while( 1 )
    {
        (*niters)++;
#if DEBUG
        printf("Iteration %d\n", *niters);
        print_matrix( phi_cur, n );
        sleep(1);
#endif
          omp_set_num_threads(TN);
        // Compute next (new) phi from current (old) phi
        for( j = 1; j < n-1; j++ )
{
#pragma omp parallel
{
#pragma omp for
            for( i = 1; i < n-1; i++ )
 {
                phi_next[ j*n+i ] = ( phi_cur[ (j-1)*n+i ] +
                                      phi_cur[ j*n+(i-1) ] +
                                      phi_cur[ j*n+(i+1) ] +
                                      phi_cur[ (j+1)*n+i ] ) / 4;
 }
// parallel region finished
  }
  }
        // If converged, we are done
        conv = 1;
        for( j = 1; j < n-1; j++ )
            for( i = 1; i < n-1; i++ )
                if( fabs( phi_next[ j*n + i ] - phi_cur[ j*n + i ] ) > THRESHOLD )
                    conv = 0;
        if ( conv )
            break;
        // Otherwise, swap pointers and continue
        tmp = phi_cur;
        phi_cur = phi_next;
        phi_next = tmp;
    }

    free( phi_cur );
    free( phi_next );
}



double compute_single_region( int n, int *niters )
{
    double *phi_cur, *phi_next, *tmp;
    int conv;
    int i, j;

    phi_cur  = (double *) malloc ( n * n * sizeof(double) );
    phi_next = (double *) malloc ( n * n * sizeof(double) );

    init_phi(phi_cur, n);
    init_phi(phi_next, n);
    // initialization of convergence criteria

    *niters = 0;
    while( 1 )
    {
        (*niters)++;
#if DEBUG
        printf("Iteration %d\n", *niters);
        print_matrix( phi_cur, n );
        sleep(1);
#endif
        conv = 1;
        // Compute next (new) phi from current (old) phi

        #pragma omp parallel reduction(min:conv)
{
#pragma omp for private(i)
        for( j = 1; j < n-1; j++ )
            for( i = 1; i < n-1; i++ )
                phi_next[ j*n+i ] = ( phi_cur[ (j-1)*n+i ] +
                                      phi_cur[ j*n+(i-1) ] +
                                      phi_cur[ j*n+(i+1) ] +
                                      phi_cur[ (j+1)*n+i ] ) / 4;

        // If converged, we are done

#pragma omp for private(i)
        for( j = 1; j < n-1; j++ )
            for( i = 1; i < n-1; i++ )
                if( fabs( phi_next[ j*n + i ] - phi_cur[ j*n + i ] ) > THRESHOLD )
                    conv = 0;
}
        if ( conv )
            break;
        // Otherwise, swap pointers and continue
        tmp = phi_cur;
        phi_cur = phi_next;
        phi_next = tmp;
    }

    free( phi_cur );
    free( phi_next );
	return *phi_next;
}




int main( int argc, char *argv[] )
{
    int n, niters;
    double t0, t1, ts, tt, res;
    double speed=ts/tt;
    int A=TN;
    if (argc != 2)
    {
        fprintf( stderr, "Usage: %s <n>\n", argv[0] );
        exit( -1 );
    }

    n = atoi( argv[1] );

    printf("Version             | #iters | #threads | #Time (s) | Speedup \n");
    printf("------------------------------------------------------------- \n");

    t0 = omp_get_wtime();
// the reason why we input nilters's address is to use it as global and local variable.
    res=compute_seq( n, &niters );
    t1 = omp_get_wtime();
    ts = t1 - t0;
    printf("Sequential          | %6d | %8d | %9.2f | %7.2f \n", niters, 1, ts, 1.0);
//printf("\n\n phi=%e\n\n",res);
    // Add your code here!
// outer
    t0 = omp_get_wtime();
    compute_outer( n, &niters );
    t1 = omp_get_wtime();
    tt = t1 - t0;
speed=ts/tt;
    printf("Parallel_out        | %6d | %8d | %9.2f | %7.2f \n", niters, A, tt, speed);

// inner
    t0 = omp_get_wtime();
    compute_inner( n, &niters );
    t1 = omp_get_wtime();
    tt = t1 - t0;
speed=ts/tt;
    printf("Parallel_in         | %6d | %8d | %9.2f | %7.2f \n", niters, A, tt, speed);

// single_region
omp_set_num_threads(TN);
    t0 = omp_get_wtime();
    res=compute_single_region( n, &niters );
    t1 = omp_get_wtime();
    tt = t1 - t0;
speed=ts/tt;
    printf("Signle_region       | %6d | %8d | %9.2f | %7.2f \n", niters, A, tt, speed);

printf("\n\n phi=%e\n\n",res);






    printf("------------------------------------------------------------- \n");

    return 0;
}

