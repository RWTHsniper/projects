
// mpicc ex1.c -Wall -o ex1.x
// mpiexec -np 4 ./ex1.x 4


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"


int main( int argc, char *argv[] )
{   
    int myrank, nprocs;
    int n, mat_size;
    int i, j, tag=0;
 //srand(time(NULL));
//srand(myrank*100);
 //   MPI_Request request;

    MPI_Init(&argc,&argv);

 if ( argc != 2 )
    {
        fprintf( stderr, "Usage: %s n\n", argv[0] );
        exit( -1 );
    }

    n = atoi(argv[1]);
    mat_size = (n/2) * (n/2);
    float sm1[mat_size][mat_size], sm2[mat_size][mat_size]; // SM: submatrix, sm1 for sending, sm2 for receiving.

    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

printf("[%d] I am printing sm1\n",myrank);
for (i=0;i<mat_size;i++)
{	for (j=0;j<mat_size;j++)
{ 
sm1[i][j]=rand() % 100*myrank+3;
printf("%f", sm1[i][j]);
}
printf("\n");
}

if (myrank==0)
printf("\n");


if (myrank == 0 || myrank ==2)
{
MPI_Ssend(sm1, mat_size*mat_size, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
MPI_Recv(sm2, mat_size*mat_size, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
for (i=0;i<mat_size;i++)
{	for (j=0;j<mat_size;j++)
{ 
sm1[i][j]=sm1[i][j]+sm2[i][j];
}
}
}
else if (myrank == 1 || myrank ==3)
{
MPI_Recv(sm2, mat_size*mat_size, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Ssend(sm1, mat_size*mat_size, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
for (i=0;i<mat_size;i++)
{	for (j=0;j<mat_size;j++)
{ 
sm1[i][j]=sm2[i][j]-sm1[i][j];
}
}
}

// Synchronization
MPI_Barrier( MPI_COMM_WORLD);

printf("[%d] I am printing Bsm1\n",myrank);
for (i=0;i<mat_size;i++)
{	for (j=0;j<mat_size;j++)
{ 
printf("%f", sm1[i][j]);
}
printf("\n");
}

if (myrank==0)
printf("\n");

// Go for C

if (myrank == 0 || myrank ==1)
{
MPI_Ssend(sm1, mat_size*mat_size, MPI_FLOAT, myrank+2, tag, MPI_COMM_WORLD);
MPI_Recv(sm2, mat_size*mat_size, MPI_FLOAT, myrank+2, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
for (i=0;i<mat_size;i++)
{	for (j=0;j<mat_size;j++)
{ 
sm1[i][j]=sm1[i][j]+sm2[i][j];
}
}
}
else if (myrank == 2 || myrank ==3)
{
MPI_Recv(sm2, mat_size*mat_size, MPI_FLOAT, myrank-2, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Ssend(sm1, mat_size*mat_size, MPI_FLOAT, myrank-2, tag, MPI_COMM_WORLD);
for (i=0;i<mat_size;i++)
{	for (j=0;j<mat_size;j++)
{ 
sm1[i][j]=sm2[i][j]-sm1[i][j];
}
}
}

// Synchronization
MPI_Barrier( MPI_COMM_WORLD);

printf("[%d] I am printing Csm1\n",myrank);
for (i=0;i<mat_size;i++)
{	for (j=0;j<mat_size;j++)
{ 
//sm1[i][j]=rand() % 50;
printf("%f", sm1[i][j]);
}
printf("\n");
}




MPI_Finalize();


return 0;

}
