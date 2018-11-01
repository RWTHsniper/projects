// DO NOT MODIFY THIS FILE
#ifndef INTERFACES_H
#define INTERFACES_H

#define NEN 4
#define NDF 4
#define NSD 3

// Parameters
extern const unsigned nen;
extern const unsigned ndf;
extern const unsigned nsd;
extern const unsigned nn_fine;
extern const unsigned nn_coarse;
extern const unsigned ne_fine;
extern const unsigned ne_coarse;

// I/O routines
void read_nodes(double* mxyz, unsigned nn, unsigned nsd, char* file_name);
void read_data(double* data, unsigned nn, unsigned ndf, char* file_name);
void read_connectivity(int* mien, unsigned ne, unsigned nen, char* file_name);
void write_data(double* data, unsigned nn, unsigned ndf, char* file_name);
void swapbytes(char* array, int nelem, int elemsize);

// Projection routines
int check_with_tolerance(double* xi_eta_zeta, const double* const xe, const double* const ye, const double* const ze,
                         const double xx, const double yy, const double zz, double tol);
double interpolate_data(double* xi_eta_zeta, double* coarse_data);

// Checks
void test_check_with_tolerance(int retValue, double* xi_eta_zeta, const double* const xe,
                               const double* const ye, const double* const ze,
                               const double xx, const double yy, const double zz, double tol);
void check_input(double (*mxyz_coarse)[NSD], int (*mien_coarse)[NEN], double (*data_coarse)[NDF], double (*mxyz_fine)[NSD], unsigned nn_fine, unsigned nn_coarse, unsigned ne_coarse);
void check_output(double (*data_fine)[NDF], unsigned nn_fine, double tol);
void read_nodes_correct(double* mxyz, unsigned nn, unsigned nsd, char* file_name);
void read_data_correct(double* data, unsigned nn, unsigned ndf, char* file_name);
void read_connectivity_correct(int* mien, unsigned ne, unsigned nen, char* file_name);
void write_data_correct(double* data, unsigned nn, unsigned ndf, char* file_name);

#endif // INTERFACES_H
