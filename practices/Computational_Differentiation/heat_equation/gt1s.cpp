#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cfloat>
using namespace std;
#include "dco.hpp"
using namespace dco;
typedef gt1s<double>::type DCO_TYPE;

#include "f.hpp"

int main(int argc, char* argv[]){
  cout.precision(15);
  if (argc!=3) {
    cerr << "2 parameters expected:" << endl
	 << "  1. number of spatial finite difference grid points" << endl
	 << "  2. number of implicit Euler steps" << endl;
    abort();
  }
  size_t n=atoi(argv[1]), m=atoi(argv[2]); 
	// diffusitivity coefficients.
  vector<double> c(n);
// initialize c and T
  for (size_t i=0;i<n;i++) c[i]=0.01;
  vector<double> T(n);
  for (size_t i=0;i<n-1;i++) T[i]=300.;
  T[n-1]=1700.;
  ifstream ifs("O.txt");
  vector<double> O(n);
  for (size_t i=0;i<n;i++) ifs >> O[i] ;
  vector<double> dvdc(n,0);
  // ca and Ta will be used for AD. It is necessary to define them using DCO_TYPE
  vector<DCO_TYPE> ca(n);
  vector<DCO_TYPE> Ta(n);
  // norm of difference between T_simulation and O
  DCO_TYPE va;

  double gradient_norm = 9999;
  int iter = 0;
  const int maxiter = 5e6;

  while (gradient_norm > 1e-5 && iter < maxiter) {
	  ++iter;

  for(size_t i=0;i<n;i++) {
// Load ca and Ta
// c and T are just vector. ca and Ta are special format.
	  for(size_t j=0;j<n;j++) { ca[j]=c[j]; Ta[j]=T[j]; }
	  derivative(ca[i])=1.0;
	  f(ca,m,Ta,O,va); // vec<dco>, size_t, vec<dco>, vec<double>, dco
	  dvdc[i]=derivative(va);
  }

	  gradient_norm = 0.0;
	  for(size_t j=0;j<n;j++) {
		gradient_norm += dvdc[j] * dvdc[j];
	  }
	  gradient_norm = sqrt(gradient_norm);

//	  const double alpha = 1e-9;
	  const double alpha = 3e-2;
	  for(size_t j=0;j<n;j++) {
		c[j] = c[j] - alpha * dvdc[j];
	  }

	  cout << "Iteration " << iter << " gradient norm " << gradient_norm << endl;

}
  cout.precision(15);

  for(size_t i=0;i<n;i++) 
    cout << "c[" << i << "]=" << c[i] << endl;

  for(size_t i=0;i<n;i++) 
    cout << "dvdc[" << i << "]=" << dvdc[i] << endl;

  return 0;
}



