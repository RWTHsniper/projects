
double coo_to_csr(int* Ic, int* I, int* Jc, int* J_csr, int* J, double* valc, double* val, int M, int nz, int symm);

void matrix_vector(int* J_csr, int* Jc, double* valc, double*Xv, double* y, int M);

double dot_p(double* a, double* b, int k);

double getkrylov(int* I, int* J, double* val, int* Jc, int* J_csr, double* valc, double* H, int* H_csc, double* V, double* V_c,  double* w, const int j, int M, int N, int nz, double* g, double* c, double* s, int symm, int* precond, int* precond_csr, double* precond_val, int precondition);

double norm(double* v, int k);

void H_inv_g(double* H, int* H_csc, int M, int H_size, double* g, double* h_inv_g);

void dense_matrix_vector(double* A, int M, int N, double* x, double* y);

int Precondition(int* Jc, int* J_csr, double* valc, int M, int nz, int* precond, int* precond_csr, double* precond_val, int precondition);



