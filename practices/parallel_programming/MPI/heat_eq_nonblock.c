
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define threshold 0.001


int main( int argc, char *argv[] )
{
    int myrank, nprocs;
    int n, nr, tag=0;
	double t1, t2;
    int i, j, root=0, flag=1, flag2=0, count=0, count2=0;
MPI_Request request;
// Input value is size of matrix
 if ( argc != 2 )
    {
        fprintf( stderr, "Usage: %s n\n", argv[0] );
        exit( -1 );
    }

    n = atoi(argv[1]);
   // mat_size = n * n;


    MPI_Init(&argc,&argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

    int ngather[nprocs];

if (myrank ==root)
t1= MPI_Wtime();

if (myrank < nprocs-1)
{
nr= (n/nprocs);
}
else
{
nr= n-(nprocs-1)*(int)(n/nprocs);
}
int mydata=n*nr;

//printf("[%d] nr= %d\n",myrank, nr);

MPI_Barrier(MPI_COMM_WORLD);

float u[nr][n], u_new[nr][n], temp[n], temp2[n];
float u_t[n*n];

// Initializing values

for (i=0;i<nr;i++)
{
for (j=1; j<n-1;j++)
{
u[i][j]=50;
}
}


// setting up boundary condition
// Bottom
if (myrank == nprocs-1)
{
for (i=0;i<n;i++)
{
u[nr-1][i]=0;
u_new[nr-1][i]=0;
}
}
// LEFT
for (i=0;i<nr;i++)
{
u[i][0]=100;
u_new[i][0]=100;
}
// Right
for (i=0;i<nr;i++)
{
u[i][n-1]=100;
u_new[i][n-1]=100;
}
// Top
if (myrank == root)
{
for (i=0;i<n;i++)
{
u[0][i]=100;
u_new[0][i]=100;
}
}


MPI_Barrier(MPI_COMM_WORLD);
// check
/*
if (myrank == 1)
{printf("[%d]I want to show you what I got\n",myrank);
for (i=0;i<nr;i++)
{
for (j=0; j<n;j++)
{
printf("%f ",u[i][j]);
}
printf("\n");
}
printf("\n");
}
*/

//
while(flag2==0)
//while(count<310)
{
// Calculate next step without communication.

for (i=1;i<nr-1;i++)
{
for (j=1;j<n-1;j++)
{
u_new[i][j]=(u[i-1][j]+u[i][j-1]+u[i+1][j]+u[i][j+1])/4;
}
}

// First communication to calculate lowest row

if (myrank >root && myrank < nprocs-1)
{
for (i=1;i<n-1;i++)
{
MPI_Isend(&u[0][i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &request);
}
for (i=1;i<n-1;i++)
{
MPI_Recv(&temp[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
else if (myrank == root)
{
for (i=1;i<n-1;i++)
{
MPI_Recv(&temp[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
else if (myrank == nprocs-1)
{
for (i=1;i<n-1;i++)
{
MPI_Isend(&u[0][i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD,&request);
}
}
// calculate next step proc 0 to nprocs-2
if (nr>1)
{
if (myrank < nprocs-1)
{
for (i=1;i<n-1;i++)
{
if (myrank == root && nr == 1)
{
continue;
}
else
{
u_new[nr-1][i]=(u[nr-1][i+1]+u[nr-1][i-1]+u[nr-2][i]+temp[i])/4;
}
}
}
}

// Second communication to calculate highest row
if (myrank >root && myrank < nprocs-1)
{
for (i=1;i<n-1;i++)
{
MPI_Isend(&u[nr-1][i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &request);
}
for (i=1;i<n-1;i++)
{
MPI_Recv(&temp2[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
else if (myrank == root)
{
for (i=1;i<n-1;i++)
{
MPI_Isend(&u[nr-1][i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &request);
}
}
else if (myrank == nprocs-1)
{
for (i=1;i<n-1;i++)
{
MPI_Recv(&temp2[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
// calculate next step
if (nr>1)
{
if (myrank > root)
{
for (i=1;i<n-1;i++)
{
if (myrank == nprocs-1 && nr ==1)
{
continue;
}
else
{
u_new[0][i]=(u[0][i+1]+u[0][i-1]+u[1][i]+temp2[i])/4;
}
}
}
}

// Next step for processors having 1 nr
if (root < myrank && myrank < nprocs-1)
{
if (nr == 1)
{
for (i=1;i<n-1;i++)
{
u_new[0][i]=(u[0][i+1]+u[0][i-1]+temp[i]+temp2[i])/4;
}
}
}

for (i=0;i<nr;i++)
{
for (j=1;j<n-1;j++)
{
if (fabs(u[i][j]-u_new[i][j]) > threshold)
{flag=0;
break;
}
else
{
flag=1;
}
}
}

// renew u

for (i=0;i<nr;i++)
{
for (j=1;j<n-1;j++)
{
u[i][j]=u_new[i][j];
//printf("[%d, %d]renewing, u[%d][%d] %f\n",myrank,count,i,j, u[i][j]);
}
}


flag2=0;
MPI_Allreduce(&flag, &flag2, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
//printf("flag2 = %d",flag2);
if (flag2 == 1)
{
// Every processor should let root know how many data are they going to send.
MPI_Gather(&mydata, 1, MPI_INT, ngather, 1, MPI_INT, root, MPI_COMM_WORLD );


if (myrank ==root)
{

//printf("I successfully gathered\n");
for (i=0;i<nprocs;i++)
{printf("%d ",ngather[i]);
}
printf("\n");
// initialize u_t
for (i=0; i< nr;i++)
{
for (j=0; j<n;j++)
{u_t[count2]=u[i][j];
count2++;
}
}


for (i=1;i<nprocs;i++)
{
for (j=0;j<ngather[i];j++)
{MPI_Recv(&u_t[count2], 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
count2++;
}
}
}
else
{
for (i=0;i<nr;i++)
{for (j=0; j<n; j++)
{MPI_Isend(&u[i][j], 1, MPI_FLOAT, root, tag, MPI_COMM_WORLD, &request);
}
}
}

MPI_Barrier(MPI_COMM_WORLD);
// Finally, showing the result
/*
if (myrank == root)
{
printf("\n Total matrix\n");
for (i=0;i<n;i++)
{
for (j=0;j<n;j++)
{
printf("%3.2f ",u_t[i*n+j]);
}
printf("\n");
}
printf("count2 = %d\n",count2);
}
*/

// ending of count
}

count++;
// end of while loop
}

if (myrank ==root)
{
t2= MPI_Wtime();
printf("Elapsed time is %f\n", t2-t1);
}
MPI_Finalize();

return 0;
}

