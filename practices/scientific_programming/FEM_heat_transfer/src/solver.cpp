
#include "solver.h"
#include "dco.hpp"
using namespace dco;
using namespace std;

typedef gt1s<double>::type DCO_TYPE;

/***************************************************************************************************
 * void femSolver::solve()
 * Driver routine for the femSolver
 ***************************************************************************************************/
// we got *stngs, *nodes, *elems
void femSolver::solve()
{
    double init_T=stngs->getInitT();
    int nn = stngs->getNn();
    int max_iter_cali=100;
    double k=stngs->getD(); // Initial guess
    double dvdk = 0;
    const double alpha = 3e-2;
	double grad_tol = 1e-6;	
	double* ref_data = new double[2211];	// 2211*3
	vector<int> Index(201,0);	// Saving the data we need only.	
	// Load reference data
	ifstream ifs("rod_data.txt");
	ofstream output[3];
	output[0].open("ka.txt");	output[0] << k << endl;
	output[1].open("Ta.txt");	
	output[2].open("Ta_x.txt");	
  for (size_t i=0;i<2211;i++){  
  	for (size_t j=0; j<3; j++) {
		  if (j==0){
			  ifs >> ref_data[i];
//	  cout << ref_data[i] << endl;
}
else	// Only taking Temperature
{
  double dummy;
	ifs >> dummy;
}
  }
//  cout << endl;
  }

    // Set initial temperature value
    for (size_t i=0;i<nn;i++){   
	// Set initial T on node which does not belong to Dirichlet BC.
	if (nodes->getBCtype(i) != 1) nodes->setT(i,init_T);
    }

// Find index to store Temperature and x-coordinate
int count =0;
for (size_t i=0; i<nn;i++){
	if ((0.5-grad_tol < nodes->getY(i)) && (0.5+grad_tol > nodes->getY(i))){
		Index[count] = i;
		output[2] << nodes->getX(i) << endl;
		output[1] << nodes->getT(i) << " ";
		count++;
	}
}
	output[1] << endl;

	cout << "Initial k = " << k << endl;
    for (size_t it=0;it<max_iter_cali;it++){
	// Run FEM solver multiple times
	// setup intial values on AD terms
	for (size_t i=0;i<nn;i++) Ta[i] = nodes->getT(i);
	double gradient_norm =0;
	ka = k; 
	derivative(ka)=1.0; // seeding
	va = 0;
	// Run explicitSolver
	explicitSolver();   // It should be modified to take Ta, ka
	for (size_t i=0;i<2211;i++){
	va = va + (Ta[i] - ref_data[i])*(Ta[i] - ref_data[i]);
	}
	dvdk = derivative(va);
	gradient_norm = sqrt(dvdk*dvdk);
	k = k - alpha * dvdk;
	cout << "[" << it << "]" " gradient norm " << gradient_norm << endl;
	cout << "va = " << va << endl;
	cout << "k = " << k << endl;
	output[0] << k << endl ;

	for (size_t i=0; i<201; i++){
		output[1] << Ta[Index[i]] << " ";
	}	// Save Ta and in each step
	output[1] << endl;

	if (gradient_norm < grad_tol){
		cout << "Gradient norm is converged!! " << endl;
		break;
	}
    }
	cout << "Initial k = " << stngs->getD() << " alpha = " << alpha << " grad_tol " << grad_tol << endl;
    cout << endl << "Total execution time is " << fixed << (float)clock()/CLOCKS_PER_SEC 
	<< " seconds." << endl;

delete[] ref_data;
for (size_t i=0;i<3;i++) output[i].close();

    return;
}


/***************************************************************************************************
 * void femSolver::explicitSolver()
 * This routine solves the equation system using explicit integration.
 ***************************************************************************************************/
void femSolver::explicitSolver()
{
    int ne =stngs->getNe();
    int nen=nenTri;
    int nn = stngs->getNn();
    int f_convergence=0;
    DCO_TYPE* RHS = new DCO_TYPE[nn];
    double dt = stngs->getDt();
    int max_iter = stngs->getNIter();
    double tol = 1e-6;

    // start iteration
    for (size_t iter=0; iter<max_iter; iter++)  
    {//[M^s]{T^(s+1)}= [M^s}]{T^s}+dt*[{F^s}+{B^s}-[K^s]{T^s}]
    // Initialize parameters
    for (size_t i=0;i<nn;i++){ // M*T = M_ii * T_i
	if (nodes->getBCtype(i) != 1) {
	    DCO_TYPE T_i = Ta[i];    double mass_i= nodes->getMass(i);
	    RHS[i] = mass_i*T_i;  
	}   // RHS = M*T is done
    }   // Ned to go for RHS = RHS + dt*[F+B-K*T]

    for (size_t ie=0;ie<ne;ie++){    
	// Assemble RHS using Element arrarys
	double* element_K=elems->getKptr(ie);   // get K in ie-th element
	double* element_F=elems->getFptr(ie);   // get elem F array
	for (size_t ien=0;ien<nen;ien++)
	{// calculate RHS: M*T+dt*[F+B-K*T]
	    // Indicial notation: RHS_i = M_ii*T_i + dt * [F_i+B_i-K_ij*T_j]
	    // currently: RHS_i = M_ii*T_i
	    // rhs = [F_i+B_i-K_ij*T_j]*dt        
	    // ien : index_i, jen : index_j, Local index in an element : Global index
	    int index_i = elems->getConn(ie, ien);
	    if (nodes->getBCtype(index_i) != 1){
		DCO_TYPE rhs = dt*element_F[ien];  // dt*(F_i+B_i)
		for (size_t jen=0;jen<nen;jen++)
		{   // -dt*K_ij*T_j
			int index_j = elems->getConn(ie, jen);  DCO_TYPE T_j = Ta[index_j];
		    rhs -= dt*ka*element_K[ien*nen+jen]*T_j;
		}
		//        cout << "rhs " << rhs << endl;
		RHS[index_i]+=rhs;        
	    }
	}
    }
    // Solve the system to calculate next Temperature M*T_(n+1) = RHS_n
    // Indicial notation: M_ii * T_i = RHS_i
    for (size_t i=0;i<nn;i++)
    {
	if (nodes->getBCtype(i) != 1) 
	{   // Calculate on row which does not belong to Dirichlet BC.
	    double mass_i = nodes->getMass(i);
	    // start solving procedure on i-th row, it is possible since mass_M is diagonal.
	    DCO_TYPE T_temp= RHS[i]/mass_i;
	    // Update new temperature
	    Ta[i] = T_temp;
	    if ((isinf(value(T_temp))) || (value(isnan(T_temp)))) {
		cout << "Inf or Nan on next T value detected " << endl;
		cout << "Suggest to decrease dt " << endl;
		exit(-1);
	    }
	}
    }
}

delete[] RHS;
return;
}
