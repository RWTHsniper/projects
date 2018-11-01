// DO NOT MODIFY THIS FILE!!!

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interfaces.h"

void check_input(double (*mxyz_coarse)[NSD], int (*mien_coarse)[NEN], double (*data_coarse)[NDF], double (*mxyz_fine)[NSD], unsigned nn_fine, unsigned nn_coarse, unsigned ne_coarse){
}

void check_output(double (*data_fine)[NDF], unsigned nn_fine, double tol){
}

void test_check_with_tolerance(int retValue, double* xi_eta_zeta, const double* const xe,
                               const double* const ye, const double* const ze,
                               const double xx, const double yy, const double zz, double tol){
}

void read_nodes_correct(double* mxyz, unsigned nn, unsigned nsd, char* file_name){
}

void read_data_correct(double* data, unsigned nn, unsigned ndf, char* file_name){
}

void read_connectivity_correct(int* mien, unsigned ne, unsigned nen, char* file_name){
}

void write_data_correct(double* data, unsigned nn, unsigned ndf, char* file_name) {
}
