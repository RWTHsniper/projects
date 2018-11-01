
// mpicc ex2.c -Wall -o ex2.x
// mpiexec -np 8 ./ex2.x

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define N 20

int main( int argc, char *argv[] )
{

int rank, nprocs, i;
int randomnumbers[N];
int root = 0;
int local_max, global_max[N];

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
srand(rank);

// generate 20 randomnumbers
for (i=0; i<N; i++)
{
randomnumbers[i] = (lrand48() % 199) - 99;
}

for (i=0; i<N; i++)
{
if (local_max < randomnumbers[i])
local_max = randomnumbers[i];
}

// reduce version
// MPI_Reduce(&local_max, &globalmax

// gather version
printf("\nThis is my rank %d and my local_max %d\n",rank, local_max);

MPI_Gather(&local_max, 1, MPI_INT, global_max, 1, MPI_INT, root, MPI_COMM_WORLD);
// check if the root got right local_max.
/*
if (rank == root)
{
for (i=0;i<nprocs;i++)
printf("\nI am the root %d and my global_max %d\n",rank, global_max[i]);
}
*/
for (i=1; i<nprocs;i++)
{
if (global_max[0]<global_max[i])
{global_max[0]=global_max[i];
}
}

MPI_Barrier(MPI_COMM_WORLD);

if (rank == root)
{
printf("\nFinally, the global_max is %d\n", global_max[0]);
}


MPI_Finalize();
return 0;

}
