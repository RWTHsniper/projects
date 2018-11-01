
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
int local_max, global_max;

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

printf("\nThis is my rank %d and my local_max %d\n",rank, local_max);

// reduce version
MPI_Reduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, root, MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD);
if (rank == root)
{
printf("\nFinally, the global_max is %d\n", global_max);
}


MPI_Finalize();
return 0;

}
