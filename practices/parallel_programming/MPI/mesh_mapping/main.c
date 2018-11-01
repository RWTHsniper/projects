#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "interfaces.h"

// divideWork: Splits "total" into npes-many "chunk"s,
//             such that all chunk sizes encountered on the different
//             PEs lie within the range [chunk,chunk+1]. Hence, the maximal
//             chunk size "maxchunk" is at most "chunk+1", or "chunk" if "total"
//             is divisible by npes.


void address_swaper(void **A, void **B)
{
void *temp = *A;
*A = *B;
*B= temp;
}

void divideWork(const int total, const int mype, const int npes, int* chunk, int* offset, int* maxchunk) {
	// CODE HERE!
	int remainder = (total % npes);

	// total work is equally divided
	if (remainder ==  0)
	{
		*chunk = (int)(total/npes);
		*maxchunk = *chunk;
	}
	else
	{
	// total work is inequally divided
		if (mype < remainder)
		{
		*maxchunk = (int)(total/npes) +1;
		*chunk = *maxchunk;
		}
		else
		{
		*chunk = (int)(total/npes);
		*maxchunk = *chunk +1;
		}
	}

	// calculate offset
	if ((mype) < remainder)
	{
		*offset = mype * (*maxchunk);
		// some processors have maximum number of chunk.
		*chunk = *maxchunk;
	}
	else
	{// if remainder =1, mype = 1, offset = maxchunk + chunk
		*offset = remainder * (*maxchunk) + (mype -remainder) * (*chunk);
	}
}

int main(int argc, char** argv) {

	// Get rank (mype) and size (npes)
	int mype, npes;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &mype );
	MPI_Comm_size( MPI_COMM_WORLD, &npes );

	int prev = (npes + mype -1) % npes; // CODE HERE: Previous node in the round-robin
	int next = (mype +1) % npes; // CODE HERE: Next node in the round-robin

  // Divide coarse elements and fine mesh nodes across PEs
  int nnc_fine, nec_coarse;
  int offset_nn_fine, offset_ne_coarse;
  int max_nnc_fine, max_nec_coarse;

  divideWork(ne_coarse, mype, npes, &nec_coarse, &offset_ne_coarse, &max_nec_coarse);
  divideWork(nn_fine, mype, npes, &nnc_fine, &offset_nn_fine, &max_nnc_fine);

  // Storage:
  double (*mxyz_coarse)[NEN][NSD] = (double(*)[NEN][NSD]) malloc(sizeof(double)*nsd*nen*max_nec_coarse);
  double (*data_coarse)[NEN][NDF] = (double(*)[NEN][NDF]) malloc(sizeof(double)*ndf*nen*max_nec_coarse);
  double (*mxyz_fine)[NSD] = (double(*)[NSD]) malloc(sizeof(double)*nsd*max_nnc_fine);
  double (*data_fine)[NDF] = (double(*)[NDF]) malloc(sizeof(double)*ndf*max_nnc_fine);
  int *node_found = (int*) malloc(sizeof(int)*max_nnc_fine);

  // Send buffers
  double (*mxyz_fine_send)[NSD] = (double(*)[NSD]) malloc(sizeof(double)*nsd*max_nnc_fine);
  double (*data_fine_send)[NDF] = (double(*)[NDF]) malloc(sizeof(double)*ndf*max_nnc_fine);
  int *node_found_send = (int*) malloc(sizeof(int)*max_nnc_fine);
  int nnc_fine_send, offset_nn_fine_send;

  // Read coarse mesh
  read_data(&mxyz_coarse[0][0][0], nec_coarse, offset_ne_coarse, nsd*nen, "Data/mxyz.coarse.local");
  read_data(&data_coarse[0][0][0], nec_coarse, offset_ne_coarse, ndf*nen, "Data/data.coarse.local");

  // Read fine mesh
  read_data(&mxyz_fine[0][0], nnc_fine, offset_nn_fine, nsd, "Data/mxyz.fine");

	// Use the following variables for timing
	double start_time = MPI_Wtime(); // CODE HERE: Start time
	MPI_Request request[5]; /* request[0]:sendbuf request[1]: recvbuf */

  // Have some indicator for the nodes which are found
  for(int i = 0; i < nnc_fine; ++i) {
    node_found[i] = 0;
  }

  // Test with increasing tolerances until all nodes are interpolated
  double tol = 0.0;
  for (int rounds=0; rounds < 2; rounds++) {
    if (mype==0) {
      printf("Running with tol=%f\n", tol);
    }

    // Round-robin loop
    for (int ipes=0; ipes <= npes-1; ipes++) {

			// Irecv communication to get new data
			if ((ipes + rounds) > 0) {

				MPI_Recv(&(nnc_fine), 1, MPI_INT, prev, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


				MPI_Recv(&(mxyz_fine[0][0]), nnc_fine*nsd, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&(data_fine[0][0]), nnc_fine*ndf, MPI_DOUBLE, prev, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&(node_found[0]), nnc_fine, MPI_INT, prev, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


				printf("[%d] MPI_Recv done, message size %d \n",mype, nnc_fine);
			}



      // Loop over all fine nodes and search for matching coarse element with given tolerance
      for(int i = 0; i < nnc_fine; ++i) {
        // Skip if node was found in a previous run
        if(node_found[i]!=0) continue;

        // Save fine coordinates
        double xx = mxyz_fine[i][0]; // The fine node x-coordinate
        double yy = mxyz_fine[i][1]; // The fine node y-coordinate
        double zz = mxyz_fine[i][2]; // The fine node z-coordinate

        // Check which coarse element contains the fine node
        for(int j = 0; j < nec_coarse; ++j){
          double xi_eta_zeta[NSD];  // The parametric coordinats (can be computed in the function)
          double xe[NEN];           // The coarse element level x-coordinates
          double ye[NEN];           // The coarse element level y-coordinates
          double ze[NEN];           // The coarse element level z-coordinates

          // Save coarse element coordinates
          xe[0] = mxyz_coarse[j][0][0];
          xe[1] = mxyz_coarse[j][1][0];
          xe[2] = mxyz_coarse[j][2][0];
          xe[3] = mxyz_coarse[j][3][0];
          ye[0] = mxyz_coarse[j][0][1];
          ye[1] = mxyz_coarse[j][1][1];
          ye[2] = mxyz_coarse[j][2][1];
          ye[3] = mxyz_coarse[j][3][1];
          ze[0] = mxyz_coarse[j][0][2];
          ze[1] = mxyz_coarse[j][1][2];
          ze[2] = mxyz_coarse[j][2][2];
          ze[3] = mxyz_coarse[j][3][2];

          // Check that node can be interpolated with given tolerance
          node_found[i] = check_with_tolerance(xi_eta_zeta, xe, ye, ze, xx, yy, zz, tol);

          // Interpolate data
          if(node_found[i]){
            int k;
            for(k = 0; k < ndf; ++k) { // Loop over degrees of freedom and interpolate
              double coarse_data[NEN];
              coarse_data[0] = data_coarse[j][0][k];
              coarse_data[1] = data_coarse[j][1][k];
              coarse_data[2] = data_coarse[j][2][k];
              coarse_data[3] = data_coarse[j][3][k];
              data_fine[i][k] = interpolate_data(xi_eta_zeta, coarse_data);
            }
            break;
          }
        }
      }


			// wait until sending (mxyz, data, node_found) is finished
			if (ipes + rounds >0)
{
				MPI_Wait(&request[0], MPI_STATUS_IGNORE); 
				MPI_Wait(&request[1], MPI_STATUS_IGNORE); 
				MPI_Wait(&request[2], MPI_STATUS_IGNORE); 
			//printf("current value of request[2] = %d",request[2]);
}
			// Isend communication
			// npes -1 + round(=1) = npes is the last step
			if (ipes + rounds < npes) {
				// calculated data stored in data_fine. Send what each processor has now

				address_swaper(&mxyz_fine_send, &mxyz_fine);
				address_swaper(&data_fine_send, &data_fine);
				address_swaper(&node_found_send, &node_found);


				MPI_Isend(&(mxyz_fine_send[0][0]), nnc_fine*nsd, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &request[0]);

				// calculated data stored in data_fine. Send what each processor has now

				MPI_Isend(&(data_fine_send[0][0]), nnc_fine*ndf, MPI_DOUBLE, next, 1, MPI_COMM_WORLD, &request[1]);

				// calculated data stored in data_fine. Send what each processor has now


				MPI_Isend(&(node_found_send[0]), nnc_fine, MPI_INT, next, 2, MPI_COMM_WORLD, &request[2]);

		// To make it sure that sending nnc, offset is complete.
			if (ipes + rounds >0)
		{
			// previous step's sending nnc_fine complete
				MPI_Wait(&request[3], MPI_STATUS_IGNORE); 
				nnc_fine_send = nnc_fine;
				MPI_Isend(&(nnc_fine_send), 1, MPI_INT, next, 3, MPI_COMM_WORLD, &request[3]);

			// sending offset_nn_fine complete
				MPI_Recv(&(offset_nn_fine), 1, MPI_INT, prev, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Wait(&request[4], MPI_STATUS_IGNORE); 
				offset_nn_fine_send = offset_nn_fine;
				MPI_Isend(&(offset_nn_fine_send), 1, MPI_INT, next, 4, MPI_COMM_WORLD, &request[4]);
		}
			else
		{
				nnc_fine_send = nnc_fine;
				offset_nn_fine_send = offset_nn_fine;
				MPI_Isend(&(nnc_fine_send), 1, MPI_INT, next, 3, MPI_COMM_WORLD, &request[3]);
				MPI_Isend(&(offset_nn_fine_send), 1, MPI_INT, next, 4, MPI_COMM_WORLD, &request[4]);
		}

			}




    }
    tol += 0.2;
  }

  int nodes_found = 0;
  for (int i =0; i < nnc_fine; ++i) {
    if (node_found[i]!=0) {
      ++nodes_found;
    }
  }
  int allnodes_found;
  // CODE HERE: allnodes_fund should contain sum over all PEs. 
	MPI_Reduce( &nodes_found, &allnodes_found, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mype==0) {
    printf("Number of nodes not found: %d\n", nn_fine-allnodes_found);
  }

  // End timing before file I/O
	double end_time = MPI_Wtime(); // CODE HERE: Write end time to this variable
	MPI_Finalize();

  // Output
  write_data(&data_fine[0][0], nnc_fine, offset_nn_fine, ndf, "Data/data.fine");

  // Deallocate memory
  free(mxyz_coarse);
  free(data_coarse);
  free(mxyz_fine);
  free(data_fine);
  free(node_found);
  free(mxyz_fine_send);
  free(data_fine_send);
  free(node_found_send);

  // Timing output
  if (mype==0) {
    double dt = end_time - start_time;
    printf("TIMING OUTPUT:\n");
    printf("--------------\n");
    printf("  Needed %e seconds to project %u nodes.\n", dt, nn_fine);
  }

  return 0;
}
