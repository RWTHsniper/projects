
// mpicc ex4.c -Wall -o ex4.x
// mpiexec -np 8 ./ex4.x

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define N 3

int myMPI_Reduce(const void *sendbuf, void *recvbuf, int count, int root, MPI_Comm comm);


int main( int argc, char *argv[] )
{

int rank, nprocs, i;
double globalmin;
double sendbuf[N], *recvbuf;
int root = 0;
recvbuf = (double*) malloc(N*sizeof(double));

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

for (i=0;i<N;i++)
{
sendbuf[i]=N*rank+i;
}

printf("I am %d and I have those values\n",rank);

for (i=0;i<N;i++)
{
printf("%e\n",sendbuf[i]);
}

myMPI_Reduce((const void *)sendbuf, (void *)recvbuf,N, root, MPI_COMM_WORLD);

if (rank==root)
{globalmin=((double *)sendbuf)[0];
printf("\nFinally, Global minimum %e\n\n",globalmin);
}
free(recvbuf);
MPI_Finalize();

return 0;
}


int myMPI_Reduce(const void *sendbuf, void *recvbuf, int count, int root, MPI_Comm comm)
{
int i, j;
int rank, nprocs;
int tag=0;
double globalmin;

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
// Root receives them
if (rank==root)
{
double *temp = (double *)malloc(N*sizeof(double));


for (i=0;i<N;i++)
{
temp[i]=((double *)sendbuf)[i];
printf("root's sendbuf[%d] = %e\n",i,((double *)sendbuf)[i]);
}

// Initializing globalmin as root's minimum value
globalmin=temp[0];
for (i=1; i<N;i++)
{
if (globalmin > temp[i])
{
globalmin = temp[i];
}
}

for (i=1;i<nprocs;i++)
{
MPI_Recv(recvbuf, N, MPI_DOUBLE,i,tag,comm,MPI_STATUS_IGNORE);
// be careful about using recvbuf
for (j=0;j<N;j++)
{
printf("recvbuf[%d] = %e\n",j,((double *)recvbuf)[j]);
if (globalmin > ((double *)recvbuf)[j])
{
globalmin = ((double *)recvbuf)[j];
}
}
}
printf("The globalmin that I want to show %e\n",globalmin);
((double *)sendbuf)[0]=globalmin;
// Find out the minimum value and store it onto recvbuf
}
// Every processor sends their values to the root 
else
{
MPI_Send(sendbuf, N, MPI_DOUBLE, root,tag, comm);
for (i=0;i<N;i++)
{
printf("[%d] sendbuf[%d] = %e\n",rank,i,((double *)sendbuf)[i]);
}
}


return 0;
}



