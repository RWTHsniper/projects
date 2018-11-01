#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "interfaces.h"
#include <time.h>
#include <omp.h>

/*

 My idea.
 This is the fastest code!!!
 
 1. Go with 0.2 tolerance firstly, save the value but do not go out.
 2. if it fits into 0th tolerance, we can save it.
 3. After that, go with only 0 tolerance loop.

 */

// It takes 47 sec to finish the program. 
// gcc -O3 -fopenmp io.c main.c definitions.c projection.c tests.c -o p2.x
// ../bin/pageri386 pinv.fine
// gs pager.ps
/* Defined on definitions.c

 // DO NOT MODIFY THIS FILE!!!
 const unsigned nen = 4;
 const unsigned ndf = 4;
 const unsigned nsd = 3;
 const unsigned nn_fine = 110618;
 const unsigned nn_coarse = 26148;
 const unsigned ne_coarse = 129968;

 */

int main(int argc, char** argv) {
	// Storage:
	// All arrays are allocated as two-dimensional arrays. Although they are dynamically
	// allocated, they can be accessed as if they would be allocated statically:
	//    int mien_coarse[ne_coarse][NEN];



	int (*mien_coarse)[NEN] = (int (*)[NEN]) malloc(
			sizeof(int) * nen * ne_coarse);
	double (*mxyz_coarse)[NSD] = (double (*)[NSD]) malloc(
			sizeof(double) * nsd * nn_coarse);
	double (*data_coarse)[NDF] = (double (*)[NDF]) malloc(
			sizeof(double) * ndf * nn_coarse);
	double (*mxyz_fine)[NSD] = (double (*)[NSD]) malloc(
			sizeof(double) * nsd * nn_fine);
	double (*data_fine)[NDF] = (double (*)[NDF]) malloc(
			sizeof(double) * ndf * nn_fine);

	// Use the following variables for timing
	double start_time; // Write start time to this variable
	double end_time;   // Write end time to this variable

	// Read coarse mesh
	read_nodes(&mxyz_coarse[0][0], nn_coarse, nsd, "Data/mxyz.coarse");
	read_connectivity(&mien_coarse[0][0], ne_coarse, nen, "Data/mien.coarse");
	read_data(&data_coarse[0][0], nn_coarse, ndf, "Data/data.coarse");

	// Read fine mesh
	read_nodes(&mxyz_fine[0][0], nn_fine, nsd, "Data/mxyz.fine");

	int i, flag_tol_s = 0, flag_tol_l = 0; /*flag for smaller and larger one */
	int m, m2;

	check_input(mxyz_coarse, mien_coarse, data_coarse, mxyz_fine, nn_fine,
			nn_coarse, ne_coarse);

	// TODO: Implement the projection of the coarse grid to the fine nodes. Use the following routine:
	// int check_with_tolerance(xi_eta_zeta, xe, ye, ze, xx, yy, zz, tol);
	// This function should return 1 if the fine node was found in the element and 0 if it was not found.
	// Use the following variables as input and output:
	double xi_eta_zeta[NSD]; // The parametric coordinats (can be computed in the function)
	double xe[NEN];           // The coarse element level x-coordinates
	double ye[NEN];           // The coarse element level y-coordinates
	double ze[NEN];           // The coarse element level z-coordinates
	double xx;                // The fine node x-coordinate
	double yy;                // The fine node y-coordinate
	double zz;                // The fine node z-coordinate
	double tol;               // A tolerance which can be used in the function
	int k, l;
	double coarse_data[NEN]; // Store the coarse element-level data values in this array // YOUR CODE STARTS HERE

	// printf("I want to see NSD %d \n",NSD); // NSD=3
	// printf("I want to see NEN %d \n",NEN); NEN=4

	start_time = omp_get_wtime();

	printf("max num of procs %d\n", omp_get_num_procs());
	omp_set_num_threads(omp_get_num_procs());
//	omp_set_num_threads(6);
	// Converting 1-base system into 0-base system.
	// m, m2: temporary value
#pragma omp parallel for private(m2)
	for (m = 0; m < ne_coarse; m++) {

		for (m2 = 0; m2 < nen; m2++)
			mien_coarse[m][m2] = mien_coarse[m][m2] - 1;
	}
#pragma omp parallel
	{
#pragma omp master
		{
			printf("Num of threads %d\n", omp_get_num_threads());
		}

	}

	// I need to find out xi_eta_zeta

	// i: fine node
	// k: data's dof
	// l: element in coarse mesh

	tol = 0;
	double tol2 = 0.12;
	//double temp_interpolation[3];

	printf("Start getting into calculation\n");

#pragma omp parallel for firstprivate(tol, tol2, flag_tol_l, flag_tol_s) private(k, l, m, xx, yy, zz, xi_eta_zeta, xe, ye, ze, coarse_data) schedule(static)
	for (i = 0; i < nn_fine; i++) {
		//printf("Working on i= %d\n",i);
		xx = mxyz_fine[i][0];
		yy = mxyz_fine[i][1];
		zz = mxyz_fine[i][2];
		{ // Loop over degrees of freedom and interpolate, ndf = 4
		  // TODO: Interpolate the data_fine[i][k]. This is the k-th degree of freedom of the i-th node.
		  // In this project the number of degrees of freedom is 1. Use for interpolation the function
		  // YOUR CODE STARTS HERE

			for (l = 0; l < ne_coarse; l++) {
				// Put element node's coordinates
				for (m = 0; m < NEN; m++) {
					xe[m] = mxyz_coarse[mien_coarse[l][m]][0];
					ye[m] = mxyz_coarse[mien_coarse[l][m]][1];
					ze[m] = mxyz_coarse[mien_coarse[l][m]][2];
				}

				if (flag_tol_l == 0) {
					if (check_with_tolerance(xi_eta_zeta, xe, ye, ze, xx, yy,
							zz, tol2) == 1) {
						for (k = 0; k < ndf; k++) {

							for (m = 0; m < NEN; m++) {
								coarse_data[m] =
										data_coarse[mien_coarse[l][m]][k];
							}
							data_fine[i][k] = interpolate_data(xi_eta_zeta,
									coarse_data);

						}
						flag_tol_l = 1;
					}
				} // We have interpolation with larger tolerance but need to go for more iteration with 0 tolerance.

				if (flag_tol_l == 1) // We got that tolerance 0.2 is accepted in a specific case. Still, we neeed to narrow it down with 0 tolerance.
						{
					if (check_with_tolerance(xi_eta_zeta, xe, ye, ze, xx, yy,
							zz, tol) == 1) {
						for (k = 0; k < ndf; k++) {
							for (m = 0; m < NEN; m++) {
								coarse_data[m] =
										data_coarse[mien_coarse[l][m]][k];
							}
							data_fine[i][k] = interpolate_data(xi_eta_zeta,
									coarse_data);
						}
						flag_tol_s = 1;
						break; // breaking from l loop, since we got interpolation with 0 toleration.
					}
				}
			} // end l

			if (flag_tol_l + flag_tol_s == 0) {
				printf("Impossible to interpolate points");
				exit(-1);
			}

			flag_tol_s = 0;
			flag_tol_l = 0;

		} // end k
	} // end i
	  // end of parallel loop

	end_time = omp_get_wtime();

	// YOUR CODE STARTS HERE

	// DO NOT CHANGE THE FOLLOWING LINES

	// Output
	write_data(&data_fine[0][0], nn_fine, ndf, "Data/data.fine");
	check_output(data_fine, nn_fine, tol);

	// Deallocate memory
	free(mien_coarse);
	free(mxyz_coarse);
	free(mxyz_fine);
	free(data_coarse);
	free(data_fine);

	// Timing output
	double dt = (end_time - start_time);
	unsigned nn = nn_fine;
	printf("TIMING OUTPUT:\n");
	printf("--------------\n");
	printf("  Needed %e seconds to project %u nodes.\n", dt, nn);
	printf("The tolerance for interpolating all nodes was: %e\n", tol);

	return 0;
}

