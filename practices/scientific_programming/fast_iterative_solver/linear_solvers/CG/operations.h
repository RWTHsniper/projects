
double coo_to_csr(int* Ic, int* I, int* Jc, int* J_csr, int* J, double* valc, int *I_csc, double* val, int M, int nz, int symm);
void matrix_vector(int* J_csr, int* Jc, int* I_csc, int* I, double* val, double* valc, double* Xv, int symm, double* y, int M);
double norm(double* v, int k);
double dot_p(double* a, double* b, int k);
void dense_matrix_vector(double* A, int M, int N, double* x, double* y);
void desne_inverse_matrix_vector(double* A, int M, double* x, double*y);






