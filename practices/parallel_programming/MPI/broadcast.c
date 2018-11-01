
// mpicc ex3.c -Wall -o ex3.x
// mpiexec -np 8 ./ex3.x

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define N 3



int mybroadcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
int myrank, nprocs;
int tag=0;
int i;

MPI_Comm_rank(comm, &myrank);
MPI_Comm_size(comm, &nprocs);

if (myrank==root)
{
for (i=1;i<nprocs;i++)
{
MPI_Send(buffer,count, datatype, i, tag, comm);
}
}
else
{
MPI_Recv(buffer, count, datatype, root, tag, comm, MPI_STATUS_IGNORE);
}

return 0;
}


int main( int argc, char *argv[] )
{

int rank, nprocs, i;
int randomnumbers[N];
int root = 0;

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
srand(rank);


if (rank==root)
{
for (i=0; i<N; i++)
{
randomnumbers[i] = (lrand48() % 199) - 99;
}
printf("I will send those values\n");
for (i=0; i<sizeof(randomnumbers)/sizeof(randomnumbers[0]); i++)
{
printf("%d\n",randomnumbers[i]);
}
printf("\n");
}

mybroadcast((void *)randomnumbers, N, MPI_INT, root, MPI_COMM_WORLD );

for (i=0; i<sizeof(randomnumbers)/sizeof(randomnumbers[0]); i++)
{
printf("My rank is %d and I got it %d\n",rank,randomnumbers[i]);
}
printf("\n");
MPI_Finalize();
return 0;

}









