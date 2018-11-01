/*
 * FIsimul.cpp
 *
 *  Created on: Nov 2, 2017
 *      Author: jaeyong
 */

#include "FIsimul.h"

#define PI 3.14159265
#define debug 0
using namespace std;
using namespace Eigen;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat; // declares a Row-major sparse matrix

FI_simul::FI_simul(grid_2d& arg_grid,unsigned int arg_num_materials,unsigned int num_ela_props,unsigned int num_arg_dam_props,unsigned int* arg_control_flags) {
	// TODO Auto-generated constructor stub
	num_materials=arg_num_materials;
	num_dam_props=num_arg_dam_props;
	restartflag=0;
	//	grid=arg_grid;

	flags.resize(5,0);
	//	Create 2D vector having elastic and plastic properties
	mat_list.resize(arg_num_materials,0);
	num_ref_props = num_ela_props;
	ref_props.resize(num_ref_props,0.0);

	if (num_ela_props>0)	helpers::create_multidim_vector(ela_prop,num_materials,num_ela_props,0.0);
	else	{
		cout << "There should be definition of elastic properties  "<<endl; abort();
	}
	if (num_arg_dam_props>0)	helpers::create_multidim_vector(dam_prop,num_materials,num_dam_props,0.0);
	//	Get some useful info from arg_grid
	for (size_t i=0;i<3;i++){
		N[i]=arg_grid.N[i];
		dx[i]=arg_grid.dx[i];
	}
	size_t intnxnynz = N[0]*N[1]*N[2];
	dim=0;
	for (size_t i=0;i<3;i++){
		if (N[i]!=1){
			dim+=1;
		}
	}
	nprocs = arg_control_flags[0];

	//	Initialize spatial dimension dependent arrays
	Pavg.resize(dim*dim,0.0e0);
	Favg.resize(dim*dim,0.0e0);

	kx.resize(N[0],0.0); ky.resize(N[1],0.0); kz.resize(N[2],0.0);
	kx_v.resize(N[0],0.0); ky_v.resize(N[1],0.0); kz_v.resize(N[2],0.0);

	//	define 2d wave vectors
	def_wave_vects();

	//	Define grid point variables to be solved and initialize them
	dam_flag.resize(N[0]*N[1]*N[2],0);
	disp.resize(N[0]*N[1]*N[2]*3,0.0e0);
	F.resize(N[0]*N[1]*N[2]*dim*dim,0.0e0);
	F_prev.resize(N[0]*N[1]*N[2]*dim*dim,0.0e0);
	P.resize(N[0]*N[1]*N[2]*dim*dim,0.0e0);
	sig.resize(N[0]*N[1]*N[2]*dim*dim,0.0e0);
	Taw_polar.resize(N[0]*N[1]*N[2]*dim*dim,0.0e0);
	F_hat.resize(N[0]*N[1]*N[2]*dim*dim,0.0e0);
	compl_tmp.resize(N[0]*N[1]*N[2],0.0e0);

	//	Decidde local or nonlocal damage model to use
	if (num_arg_dam_props > 0){
		enable_nlocal = arg_control_flags[1];
		preconditioning = arg_control_flags[2];
		if (enable_nlocal){	//	Nonlocal damage model
			n_dam.reserve(intnxnynz);
			//	Initialize dft coefs
			coef_2diff.resize(intnxnynz,0.0);
			double L[] = {N[0]*dx[0],N[1]*dx[1],N[2]*dx[2]};
			fft::dft_coeffs(&coef_2diff[0],2,N,L);
		}
		else{	//	local damage model.
			grid_dam.reserve(intnxnynz);
		}
		grid_dam_v.resize(intnxnynz,0.0);
//		grid_ndam_v.resize(intnxnynz,0.0);
		grid_ndam_v=VectorXd::Zero(intnxnynz);
	}

	//	Diagonal terms of F and F_prev sould be 1.0
	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t k=0;k<N[2];k++){
				//	Initialize dof1 things like damage
				for (size_t l=0;l<dim;l++){
					//	Diagonal terms should be 1.0
					F[i*N[1]*N[2]*dim*dim+j*N[2]*dim*dim+k*dim*dim+l*dim+l]=1.0e0;
					F_prev[i*N[1]*N[2]*dim*dim+j*N[2]*dim*dim+k*dim*dim+l*dim+l]=1.0e0;
				}
			}

}

/*
 * (forward-backward central-FD approximation of grad and div operators)
 * I tried to put Vidyasagar's wave vector as well.
 */
void FI_simul::def_wave_vects(){

	double fac[3];
	double L[3] = {N[0]*dx[0],N[1]*dx[1],N[2]*dx[2]};
	for (size_t i=0;i<dim;i++) {
		fac[i] = double(PI)*2.0/(L[i]);
	}

	if (N[0]>1){
		for (size_t i=0;i<=N[0]/2;i++)	kx[i]=fac[0]*static_cast<double>(i);
		for (size_t i=N[0]/2+1;i<N[0];i++)	kx[i]=fac[0]*(static_cast<double>(i)-static_cast<double>(N[0]));
		for (size_t i=0;i<=N[0]/2;i++)	kx_v[i]= sin(fac[0]*static_cast<double>(i)*dx[0])/dx[0];
		for (size_t i=N[0]/2+1;i<N[0];i++)	kx_v[i] = sin(fac[0]*(static_cast<double>(i)-static_cast<double>(N[0]))*dx[0])/dx[0];
	}
	else{
		kx[0] = 0.0; kx_v[0]=0.0;
	}
	if (N[1]>1){
		for (size_t i=0;i<=N[1]/2;i++)	ky[i]=fac[1]*static_cast<double>(i);
		for (size_t i=N[1]/2+1;i<N[1];i++)	ky[i]=fac[1]*(static_cast<double>(i)-static_cast<double>(N[1]));
		for (size_t i=0;i<=N[1]/2;i++)	ky_v[i]= sin(fac[1]*static_cast<double>(i)*dx[1])/dx[1];
		for (size_t i=N[1]/2+1;i<N[1];i++)	ky_v[i] = sin(fac[1]*(static_cast<double>(i)-static_cast<double>(N[1]))*dx[1])/dx[1];
	}
	else{
		ky[0]=0.0; ky_v[0]=0.0;
	}
	if (N[2]>1){
		for (size_t i=0;i<=N[2]/2;i++)	kz[i]=fac[2]*static_cast<double>(i);
		for (size_t i=N[2]/2+1;i<N[2];i++)	kz[i]=fac[2]*(static_cast<double>(i)-static_cast<double>(N[2]));
		for (size_t i=0;i<=N[2]/2;i++)	kz_v[i]= sin(fac[2]*static_cast<double>(i)*dx[2])/dx[2];
		for (size_t i=N[2]/2+1;i<N[2];i++)	kz_v[i] = sin(fac[2]*(static_cast<double>(i)-static_cast<double>(N[2]))*dx[2])/dx[2];
	}
	else{
		kz[0]=0.0; kz_v[0]=0.0;
	}

	/*	debugging
	 *
	 * 	cout << "Printing wave vectors for a check" << endl;
	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++){
			cout <<"i, j "<<i<<" "<<j<<endl;
			cout <<"kx "<< kx[i]<<" "<<ky[j]<<endl;
			cout <<"kx_v "<< kx_v[i]<<" "<<ky_v[j]<<endl;
		}
	 *
	 */

	//	1st part
	//	2nd part
	//	Kabel's notation
	//	Kabel's definition
	/*
	for (size_t i=0;i<N[0];i++)	kx[i]= 2.0e0*static_cast<double>(PI)*static_cast<double>(i)/L[0];
	for (size_t i=0;i<N[1];i++)	ky[i]= 2.0e0*static_cast<double>(PI)*static_cast<double>(i)/L[1];
	 */

	//	Let's keep working on it.
}


void FI_simul::def_accu_green(){
	//	green_h.resize(81,0.0);	// 3D case I guess.
	//	green_h.resize(N[0]*N[1]*16,0.0);	// 2D case I guess.
	//	green_h.resize(N[0]*N[1]*N[2]*pow(dim,4),0.0);
	inv_accu_h.resize(N[0]*N[1]*N[2]*dim*dim,0.0);
	//	accu_h.resize(N[0]*N[1]*N[2]*dim*dim,0.0);
	// Firstly I start with defining inverse of Accustic tensor
	double lam0 =ref_props[1]; double mu0 = ref_props[0];
	//	Go through every grid point (i,j)
	size_t N2[]={N[0],N[1],N[2],dim,dim,dim,dim};
	helpers::indexing in5(N2,5);	// Only pass first 5
	helpers::indexing in7(N2,7);	// Pass entire N2
	//	helpers::indexing in5(&N2[0],5);	// Only pass first 5
	//	helpers::indexing in7(&N2[0],7);	// Pass entire N2
	cout<<"Initialization of accu+green done"<<endl;
	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t kk=0;kk<N[2];kk++){
				//	Loop over entire domain
				double xi2[] = {kx[i],ky[j],kz[kk]}; double xi2_2 = xi2[0]*xi2[0]+xi2[1]*xi2[1]+xi2[2]*xi2[2];
				double tmpin[dim][dim], tmpo[dim][dim];	// 2D case
				if (xi2_2 > 0.0e0){// Only do it for where there is nonzero wave vector
					//	Building up inverse of accustic tensor
					for (size_t k=0;k<dim;k++)
						for (size_t l=0;l<dim;l++){
							//		2nd C0 that M.Kabel suggested
							tmpin[k][l]=lam0*xi2[k]*xi2[l]+2.0e0*mu0*xi2_2*double(k==l);
							//		First C0 that Julian uses
							//					tmpin[k][l]=lam0*xi2[k]*xi2[l]+mu0*(xi2_2*double(k==l)+xi2[k]*xi2[l]);
						}
					/*
				cout << "wave vector "<< xi2[0] <<" "<<xi2[1]<<endl;
				cout << "lam0 and mu0 "<<lam0<<" "<<mu0<<endl;
					 */
					//	Make up accustic tensor
					math_op::Inv_M(&tmpo[0][0],&tmpin[0][0],dim);
					/*
					 * Debugging...
				cout << "building up Accustic tensor. I, J KK "<< i<<" "<<j<<" "<<kk<<endl;
				helpers::print_matrix(&tmpin[0][0],dim,dim);
				cout <<"print Inv_M " <<endl;
				helpers::print_matrix(&tmpo[0][0],dim,dim);
					 */
					for (size_t k=0;k<dim;k++)
						for (size_t l=0;l<dim;l++){
							//						cout << "index "<<i<<" "<<j<<" "<<kk<<" "<<k<<" "<<l<<" "<<in5.g(i,j,kk,k,l)<< endl;
							inv_accu_h[in5.g(i,j,kk,k,l)]=tmpo[k][l];
						}
					//cout << "start making green operator"<<endl;
					//	Make up	Green operator
					/*
				for (size_t k=0;k<dim;k++)
					for (size_t L=0;L<dim;L++)
						for (size_t m=0;m<dim;m++)
							for (size_t NN=0;NN<dim;NN++){
								double km = double(k==m);
								double term1 = km*xi2[L]*xi2[NN]/(2.0e0*mu0*xi2_2);
								double term2 = -lam0*xi2[k]*xi2[L]*xi2[m]*xi2[NN]/(2.0e0*mu0*(lam0+2.0e0*mu0)*xi2_2*xi2_2);
				green_h[in7.g(i,j,kk,k,L,m,NN)] = term1 + term2;
//	Let's try different method to calculate green_h
								green_h[in7.g(i,j,kk,k,L,m,NN)] =xi2[L]*xi2[NN]*accu_h[in5.g(i,j,kk,k,m)];
								}
					 */

					//				cout << "finish making green operator"<<endl;
				}//		end	if ((xi2[0]!=0)&&(xi2[1]!=0)){// Only do it for where there is nonzero wave vector
			}
	cout<<"Finishing accu+green done"<<endl;
}

void FI_simul::save_prev_F(){

	size_t len_arr =N[0]*N[1]*N[2]*dim*dim;
	memcpy(&F_prev[0],&F[0],sizeof(F[0])*len_arr);
}

/*
 * Initialize Piola-Kirchhoff stress field according to given initial F field.
 */
void FI_simul::calc_P_field(grid_2d& grid){

	//	cout << "getting into initialization of Piola-Kirchhoff tensor"<<endl;

	size_t N1[]={N[0],N[1],N[2],dim,dim,dim,dim};
	helpers::indexing in3(N1,3);	// Only pass first 3
	helpers::indexing in5(N1,5);	// Only pass first 5
	unsigned int* arr_mat_indx=grid.get_mat_indx();
	//	Temporary array for 2D case
	double FF[3][3],PP[3][3];	// Used for calculating stress
	FF[0][2]=0.0;FF[1][2]=0.0;FF[2][2]=1.0;FF[2][0]=0.0;FF[2][1]=0.0;


	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t kk=0;kk<N[2];kk++){
				//	Take material params
				unsigned int mat_indx = arr_mat_indx[in3.g(i,j,kk)];
				double mu=ela_prop[mat_indx][0]; double lam = ela_prop[mat_indx][1];

				if (dim==2){
					for (size_t a=0;a<dim;a++)
						for (size_t b=0;b<dim;b++){
							FF[a][b]=F[in5.g(i,j,kk,a,b)];
							PP[a][b]=P[in5.g(i,j,kk,a,b)];
						}
					//	Calculate first Piola-Kirchhoff stress.
					elasticity::hyper::neo_hook::piola_kirch_stress(&PP[0][0],&FF[0][0],mu,lam,1);
					for (size_t a=0;a<dim;a++)
						for (size_t b=0;b<dim;b++){
							P[in5.g(i,j,kk,a,b)]=PP[a][b];
						}
				}
				else if (dim==3){
					//	Calculate first Piola-Kirchhoff stress.
					//				cout << "Print Fgrad on each grid P"<<endl;
					//				helpers::print_matrix(&F[in5.g(i,j,kk,0,0)],3,3);
					elasticity::hyper::neo_hook::piola_kirch_stress(&P[in5.g(i,j,kk,0,0)],&F[in5.g(i,j,kk,0,0)],mu,lam,1);
				}
				//	If it is elastic-damage simulation, damage factor should be applied.
				if (num_dam_props>0){
					double fD=0.0;
					if (enable_nlocal){
						fD = n_dam[in3.g(i,j,kk)].get_fD();
					}
					else{
						fD = grid_dam[in3.g(i,j,kk)].get_fD();
					}
					for (size_t a=0;a<dim;a++)
						for (size_t b=0;b<dim;b++){
							P[in5.g(i,j,kk,a,b)] *= fD;
						}
				}
			}
	//	cout << "Finished into initialization of Piola-Kirchhoff tensor"<<endl;

}

/*	My user-defined material
 *
 */
void FI_simul::umat(grid_2d& grid){

	if (enable_nlocal) update_nloc_damage(2600, 1.0e-6);
	else update_damage();
	calc_P_field(grid);

}

/*	update_damage
 * update damage value on each grid point
 */
void FI_simul::update_damage(){

	size_t N1[]={N[0],N[1],N[2],dim,dim};
	helpers::indexing in3(N1,3);	// Only pass first 3
	helpers::indexing in5(N1,5);	// Only pass first 5

	double F33[3][3];
	math_op::initiali_arr(&F33[0][0],0.0e0,9);
	F33[2][2] = 1.0;

	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t kk=0;kk<N[2];kk++){
				size_t indx3 = in3.g(i,j,kk);
				if (dim==2){
					for (size_t k=0;k<dim;k++)
						for (size_t l=0;l<dim;l++)	{
							F33[k][l]= F[in5.g(i,j,kk,k,l)];
						}
					grid_dam[indx3].updateF(&F33[0][0]);
				}
				else{
					grid_dam[indx3].updateF(&F[in5.g(i,j,kk,0,0)]);
				}
				//				grid_dam[indx3].update_dam();	//	update damage variable
				int flag =grid_dam[indx3].update_dam();
				if	((flags[0]==0)&&(flag==1)){
					cout << "damage initiated" << endl;
					flags[0] = 1;
				}
			}
}

/*	update_nloc_damage
 *	calculate nonlocal damage effect using FFT-based scheme
 *	Remind!!!
 *
 *	compl_tmp	:	FFT(dam)
 *	Taw_polar[Nx*Ny*Nz*dim*dim] is used for saving arrays to save memory...???.
 *	Taw_polar[0:Nx*Ny*Nz-1]: dam
 *	Taw_polar[Nx*Ny*Nz:2*Nx*Ny*Nz-1]: ndam
 *	Taw_polar[2*Nx*Ny*Nz:3*Nx*Ny*Nz-1]: ri. FFT(nloc_dam) -> residual of nonlocal damage field
 *	Taw_polar[2*Nx*Ny*Nz:3*Nx*Ny*Nz-1]: ri. FFT(nloc_dam) -> residual of nonlocal damage field
 *
 */
void FI_simul::update_nloc_damage(size_t maxiter,double rel_tol){

	//	Firstly, we need to know some terms from local damage model.
	//	kloc, phid, dY_dndam(H)

	size_t N1[]={N[0],N[1],N[2],dim,dim};
	size_t inxnynz =N[0]*N[1]*N[2];
	helpers::indexing in3(N1,3);	// Only pass first 3
	helpers::indexing in5(N1,5);	// Only pass first 5
//	vector<double>	incr_nloc(inxnynz,0.0);
	VectorXd	incr_nloc=VectorXd::Zero(inxnynz);
	Eigen::setNbThreads(nprocs);


	if (preconditioning==1) cgD.setTolerance(rel_tol);
	else if (preconditioning==2) cgLU.setTolerance(rel_tol);
	else cgI.setTolerance(rel_tol);

	double F33[3][3];
	math_op::initiali_arr(&F33[0][0],0.0e0,9);
	F33[2][2] = 1.0;

	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t kk=0;kk<N[2];kk++){
				size_t indx3 = in3.g(i,j,kk);
				if (dim==2){
					for (size_t k=0;k<dim;k++)
						for (size_t l=0;l<dim;l++)	{
							F33[k][l]= F[in5.g(i,j,kk,k,l)];
						}
					n_dam[indx3].updateF(&F33[0][0]);
				}
				else{
					n_dam[indx3].updateF(&F[in5.g(i,j,kk,0,0)]);
				}
			}
	std::fill(dam_flag.begin(), dam_flag.end(), 0);
	size_t iter,ndamage=0;
	for (iter=0;iter<maxiter;iter++){
		ndamage = 0; 	//number of grid points going through damage development
		/*	step 1:
		 * set F on each grid point and build up vectors of local damage term.
		 */
		if(debug)cout<<"step1. Build up dA12,dA21,dA22,r_loc"<<endl;
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					size_t indx3 = in3.g(i,j,kk);
					//				n_dam[indx3].update_dam();	//	update damage variable
					int damflag =n_dam[indx3].update_dam();
					if(debug) cout<<"1.1. Damage and nonlocal damage just after calculating dA12...."<<endl;
					if(debug)	cout << n_dam[indx3].dam << " "<< n_dam[indx3].ndam<<endl;
					if (damflag){
						//	Check if it is very first time for damage development.
						if ((flags[0]==0)){
							cout << "damage initiated" << endl;
							flags[0] = 1;
						}
						ndamage++; dam_flag[indx3] = 1;
						dA12[indx3] = -1.0;
						dA21[indx3] = n_dam[indx3].H;
						dA22[indx3] = n_dam[indx3].kloc;
						r_phid[indx3] = n_dam[indx3].get_phid();
						//						if(debug) cout<<"da12,da21,da22,r_phid  "<<dA12[indx3]<<" "<<dA21[indx3]<<" "<<dA22[indx3]<<" "<<r_phid[indx3]<<endl;
						//	LHS		dY/dDb*(1/kloc)
					}
					else{
						dA12[indx3] = 0.0;
						dA21[indx3] = 0.0;
						dA22[indx3] = 1.0;
						r_phid[indx3] = 0.0;
					}
				}

		if (ndamage==0){
//			cout << "nloc solving is done."<<endl;
			break;
		}
		if (debug)cout<<"print r_phid"<<endl;
		if (debug)helpers::print_matrix(&r_phid[0],inxnynz,1);
		if (debug)cout<<"print dA12"<<endl;
		if (debug)helpers::print_matrix(&dA12[0],inxnynz,1);
		if (debug)cout<<"print dA21"<<endl;
		if (debug)helpers::print_matrix(&dA21[0],inxnynz,1);
		if (debug)cout<<"print dA22"<<endl;
		if (debug)helpers::print_matrix(&dA22[0],inxnynz,1);

		/*	step 2
		 * Try to solve r_nloc using FD scheme.
		 * Maybe, FFT-based residual calculation is somehow wrong.
		 */

		//	step 2.1. Calculate r_nloc = fdop[0]*Dbar - D
		if (debug)cout<< "show nloc damgage"<<endl;
		if (debug) helpers::print_matrix(&grid_ndam_v[0],inxnynz,1);
		if (debug)cout<< "show loc damgage"<<endl;
		if (debug) helpers::print_matrix(&grid_dam_v[0],inxnynz,1);
//		r_nloc=eigfdop[0].mat_vect(&r_nloc[0],&grid_ndam_v[0]);
		r_nloc=eigfdop[0]*grid_ndam_v;
		for (size_t i=0;i<inxnynz;i++) r_nloc[i] -= grid_dam_v[i];
		//		for (size_t i=0;i<inxnynz;i++) r_nloc[i] = 0.0;
		if (debug)cout<<"print r_nloc"<<endl;
		if (debug)helpers::print_matrix(&r_nloc[0],inxnynz,1);

		/*	step 2
		 * It seems weired...
		 * calculate residual of nonlocal damage field.
		 * It is done by using FFT-based differentiation
		if (debug)cout<< "start step2. Calculate residual of nonlocal damage field using FFT"<<endl;
		for (size_t i=0;i<inxnynz;i++){
			Taw_polar[i] = static_cast<complex<double> >(n_dam[i].dam);
			Taw_polar[i+inxnynz] = static_cast<complex<double> >(n_dam[i].ndam);
			if (debug)cout<<n_dam[i].dam<<" "<<n_dam[i].ndam<<endl;
		}
		//	2.1	take FFT of damage and nloc_damage field.
		if (debug)cout<< "start step2.1. do FFT of local and nonlocal damage field"<<endl;
		if (debug)cout<<"alpha "<<(dam_prop[0][4])<<endl;
		fft::fftw(N,dim,&Taw_polar[0],&Taw_polar[0],1,nprocs);
		fft::fftw(N,dim,&Taw_polar[inxnynz],&Taw_polar[inxnynz],1,nprocs);
		if (debug)cout<<"showing dam in frequency space."<<endl;
		if	(debug)helpers::print_matrix(&Taw_polar[0],N[0],N[1]);
		if (debug)cout<<"showing ndam in frequency space."<<endl;
		if	(debug)helpers::print_matrix(&Taw_polar[inxnynz],N[0],N[1]);
				//	2.2	go through for loop over every grid point to calculate residual.
		if (debug)cout<< "start step2.2. operation on every G.P."<<endl;
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					size_t indx3 = in3.g(i,j,kk);
//	At this point!!! Only consider homogeneous distribution of alph.
// Maybe someday, inhomogeneous distribution of alpha can be considered.
//
					Taw_polar[2*inxnynz+indx3]= Taw_polar[inxnynz+indx3] - Taw_polar[indx3];
					Taw_polar[indx3] = Taw_polar[inxnynz+indx3]*(static_cast<complex<double> >(1.0)-static_cast<complex<double> >(dam_prop[0][4])*coef_2diff[indx3])-Taw_polar[indx3];
				}
		//	2.3	try ifft
		if (debug)cout<< "start step 2.3. Try ifft to calculate r_nloc"<<endl;
		fft::fftw(N,dim,&Taw_polar[0],&Taw_polar[0],0,nprocs);
		fft::fftw(N,dim,&Taw_polar[2*inxnynz],&Taw_polar[2*inxnynz],0,nprocs);
		if (debug)	cout <<"For a check, printing Dbar-D"<<endl;
		if (debug) helpers::print_matrix(&Taw_polar[2*inxnynz],1,inxnynz);
		//	2.4	save values back to r_nloc
		if (debug)cout<< "start step 2.4. Retrieve r_nloc from ifft result. printing r_nloc"<<endl;
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					size_t indx3 = in3.g(i,j,kk);
					r_nloc[indx3] =Taw_polar[indx3].real();
				}
		if (debug)helpers::print_matrix(&r_nloc[0],N[0],N[1]);
		 */



		/*	step 3
		 * Schur complement method on the problem.
		 * The problem is...
		 * [fdop, dA12]*[incr_nloc] = [-r_nloc]
		 * [dA21, dA22]*[incr_loc]  = [-r_phid]
		 */
		if (debug)cout<< "start step3. Schur complement start."<<endl;
		VectorXd	diag=VectorXd::Zero(inxnynz);
		VectorXd	rs1=VectorXd::Zero(inxnynz);
		if (debug)cout<< "start step3.1. diag and rs1 are printed out."<<endl;
		for (size_t i=0;i<inxnynz;i++){
			diag[i] = dA12[i]/dA22[i]*dA21[i];
			rs1[i] = -(r_nloc[i])-dA12[i]/dA22[i]*(-r_phid[i]);
			if(debug)cout<<diag[i]<<" "<<rs1[i]<<endl;
		}
		if (debug)cout<< "start step3.1.2. Operation matrix before adding diagonal terms."<<endl;
		//		if (debug)fdop[0].print_csr();
//		fdop[0].add_diag(&diag[0],1.0,-1.0);
		for(size_t i=0;i<inxnynz;i++) eigfdop[0].coeffRef(i,i) -= diag[i];
		if (debug)cout<< "start step3.2. Solve CG. showing operator matrix."<<endl;
		if (preconditioning==1){ 
			//		for (size_t i=0;i<inxnynz;i++) incr_nloc[i] *= 0.5;
//			fdop[0].PCG(&incr_nloc[0],&rs1[0],rel_tol);	//	Increment of nonlocal damage is solved.
//			incr_nloc = cgD.compute(eigfdop[0]).solve(rs1);
			incr_nloc = cgD.compute(eigfdop[0]).solveWithGuess(rs1,incr_nloc);
		}
		else if (preconditioning==2){
			incr_nloc = cgLU.compute(eigfdop[0]).solveWithGuess(rs1,incr_nloc);
		}
		else	incr_nloc = cgI.compute(eigfdop[0]).solve(rs1);	//	Increment of nonlocal damage is solved.
		//		if (debug)fdop[0].print_csr();
//		fdop[0].set_values(&savefdop[0]);	//	reload original values of fdop matrix.
		for(size_t i=0;i<inxnynz;i++) eigfdop[0].coeffRef(i,i) = savefdop[i];	//	reload original values of fdop matrix.
		if (debug)cout<< "start step3.3. Showing increment of nonlocal damage."<<endl;
		if (debug) helpers::print_matrix(&incr_nloc[0],inxnynz,1);
		if (debug)cout<< "start step4. Save local and nonlocal damage."<<endl;
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					size_t indx3 = in3.g(i,j,kk);
					//					grid_ndam_v[indx3] += incr_nloc[indx3];
					n_dam[indx3].ndam += incr_nloc[indx3];
					double dlam = 1.0/dA22[indx3]*(-(r_phid[indx3])-dA21[indx3]*incr_nloc[indx3]);
					if(dlam<-1.0e-4){
						cout <<"Something is wrong on dlam. dlam = "<<dlam<<endl;
						cout <<"r_phi, dA21, incr_nloc "<<r_phid[indx3]<<" "<<dA21[indx3]<<" "<<incr_nloc[indx3]<<dlam;
						//						abort();
					}
					n_dam[indx3].dam += dlam;
					n_dam[indx3].xid += dlam;
					grid_dam_v[indx3] =n_dam[indx3].dam;
					grid_ndam_v[indx3] =n_dam[indx3].ndam;
					if ((abs(n_dam[indx3].dam)>1.0)&&(abs(n_dam[indx3].xid)>1.0)&&(abs(n_dam[indx3].ndam)>1.0)){
						cout << "Damage or xid or nonlocal damage exceeding 1.0" << endl;
						cout << "i,j,k "<<i<<" "<<j<<" "<<kk<<endl;
						cout <<"damage, xid, ndam "<<n_dam[indx3].dam<<" "<<n_dam[indx3].xid<<" "<<n_dam[indx3].ndam<<endl;
						abort();
					}
					if (((n_dam[indx3].dam)<0.0)&&((n_dam[indx3].xid)<0.0)&&(abs(n_dam[indx3].ndam)<0.0)){
						cout << "Damage or xid or nonlocal damage less than 0.0" << endl;
						cout << "i,j,k "<<i<<" "<<j<<" "<<kk<<endl;
						cout <<"damage, xid, ndam "<<n_dam[indx3].dam<<" "<<n_dam[indx3].xid<<" "<<n_dam[indx3].ndam<<endl;
						abort();
					}
					//	Increment of damage	      dv2(j) = 1.0d0/dA22(j)*(A2(j)-dA21(j)*dv1(j))
					//					grid_dam_v[i] += 1.0/dA22[i]*(-(r_phid[i])-dA21[i]*incr_nloc[i]);
					if (debug)cout<<"local and nonlocal and xid  "<<n_dam[indx3].dam<<" "<<n_dam[indx3].ndam<<" "<<n_dam[indx3].xid<<endl;
				}
	}
	if ((iter == maxiter)&&(ndamage!=0)){
		cout << "Nonlocal damage solver could not converge. Abort the program.	"<<endl;
		cout << "ndamage = 	"<<ndamage<<endl;
		abort();
	}
}

/*
 * save_damage
 * save damage variables for visualization.
 * Put nonlocal damage save using enable_nlocal flag
 */
void FI_simul::save_damage(){

	size_t N1[]={N[0],N[1],N[2]};
	helpers::indexing in3(N1,3);	// Only pass first 3

	if (enable_nlocal){
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					size_t indx3 = in3.g(i,j,kk);
					grid_dam_v[indx3] = n_dam[indx3].get_dam();	//	save damage variable
					grid_ndam_v[indx3] = n_dam[indx3].get_ndam();	//	save nonlocal damage variable
				}
	}
	else{
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					size_t indx3 = in3.g(i,j,kk);
					grid_dam_v[indx3] = grid_dam[indx3].get_dam();	//	save damage variable
				}
	}
}

void FI_simul::calc_P_F_avg(){
	//		math_op::initiali_arr(&Pavg[0][0],0.0e0,9);
	size_t N1[]={N[0],N[1],N[2],dim,dim};
	double nxnynz = double(N[0]*N[1]*N[2]);
	helpers::indexing in5(N1,5);	// Only pass first 5
	for (size_t a=0;a<dim;a++)
		for (size_t b=0;b<dim;b++){
			Pavg[a*dim+b]=0.0e0;
			Favg[a*dim+b]=0.0e0;
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++)
					for (size_t k=0;k<N[2];k++){
						Pavg[a*dim+b]+=P[in5.g(i,j,k,a,b)];
						Favg[a*dim+b]+=F[in5.g(i,j,k,a,b)];
					}
		}
	for (size_t a=0;a<dim;a++)
		for (size_t b=0;b<dim;b++){
			Pavg[a*dim+b]/= nxnynz;
			Favg[a*dim+b]/= nxnynz;
		}
}

void FI_simul::convert_P_2_sig(){
	size_t N1[]={N[0],N[1],N[2],dim,dim};
	helpers::indexing in5(N1,5);	// Only pass first 3

	if (dim==3){
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t k=0;k<N[2];k++){
					size_t index = in5.g(i,j,k,0,0);
					elasticity::convert_P_2_Sig(&sig[index],&P[index],&F[index]);
				}
	}
	else if (dim==2){
		double P_gp[3][3],F_gp[3][3],sig_gp[3][3];
		//	Initialize temp arrays
		math_op::initiali_arr(&P_gp[0][0],0.0,9);
		math_op::initiali_arr(&F_gp[0][0],0.0,9);
		math_op::initiali_arr(&sig_gp[0][0],0.0,9);
		F_gp[2][2]= 1.0;	// Prevent singular tensor.
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t k=0;k<N[2];k++){
					//	Load values onto 3*3 arrays
					for (size_t ii=0;ii<dim;ii++)
						for (size_t jj=0;jj<dim;jj++){
							size_t ind2 = in5.g(i,j,k,ii,jj);
							P_gp[ii][jj]= P[ind2];
							F_gp[ii][jj]= F[ind2];
							sig_gp[ii][jj]= sig[ind2];
						}
					//	Load every tensor onto *_gp tensor (3*3) and do transformation.
					elasticity::convert_P_2_Sig(&sig_gp[0][0],&P_gp[0][0],&F_gp[0][0]);
					for (size_t ii=0;ii<dim;ii++)
						for (size_t jj=0;jj<dim;jj++){
							sig[in5.g(i,j,k,ii,jj)] = sig_gp[ii][jj];
						}
				}


	}
}
/*
 * 	This member calculated Polarization field in frequency space
 *
 */

void FI_simul::polar_field(){

	size_t N1[]={N[0],N[1],N[2],dim,dim,dim,dim};
	helpers::indexing in3(N1,3);	// Only pass first 3
	helpers::indexing in5(N1,5);	// Only pass first 5
	helpers::indexing in7(N1,7);	// Only pass first 7
	/*	Presteps
	 * load prev_F, arr_mat_indx
	 */
	save_prev_F();	// I need this array to compare values
	//	Get mat_indx array
	/*	step 1	calculate P and Polarization field
	 *
	 */
	double mu0=ref_props[0]; double lam0 = ref_props[1];

	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t kk=0;kk<N[2];kk++){
				//	Calculate polarization field T = P - C0:F
				double C0F[dim][dim];
				for (size_t k=0;k<dim;k++)
					for (size_t l=0;l<dim;l++){
						double trF=math_op::trace(&F[in5.g(i,j,kk,0,0)],dim);
						C0F[k][l]=static_cast<double>(k==l)*(lam0*trF)+2.0e0*mu0*F[in5.g(i,j,kk,k,l)];
						Taw_polar[in5.g(i,j,kk,k,l)]=(P[in5.g(i,j,kk,k,l)]-C0F[k][l]);
					}
			}
	//	For each grid point, take FFT of Taw_polar[k][l]
	for (size_t k=0;k<dim;k++)
		for (size_t l=0;l<dim;l++){
			//	Try Taw_polar[k][l] FFT
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++)
					for (size_t kk=0;kk<N[2];kk++){
						compl_tmp[in3.g(i,j,kk)]=Taw_polar[in5.g(i,j,kk,k,l)];
					}
			//	Try FFT(Taw_polar[k][l] field)
			fft::fftw(N,dim,&compl_tmp[0],&compl_tmp[0],1,nprocs);
			//	Rewriting FFT done values onto Taw_polar[k][l]
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++)
					for (size_t kk=0;kk<N[2];kk++){
						Taw_polar[in5.g(i,j,kk,k,l)]=compl_tmp[in3.g(i,j,kk)];
					}
		}
}


/*
 * Run forward FFT and inverse FFT solve for deformation gradient field.
 * Before calling this member, P field should be updated.
 */

void FI_simul::run_simul(double* F_step){

	double nxnynz = double(N[0]*N[1]*N[2]);
	size_t N1[]={N[0],N[1],N[2],dim,dim,dim,dim};
	helpers::indexing in3(N1,3);	// Only pass first 3
	helpers::indexing in5(N1,5);	// Only pass first 5
	helpers::indexing in7(N1,7);	// Only pass first 7
	/*	Presteps
	 * load prev_F, arr_mat_indx
	 */
	save_prev_F();	// I need this array to compare values
	//	Get mat_indx array
	/*	step 1	calculate and Polarization field
	 *
	 */

	/*	step 2	Take FFT(Taw_polar)
	 * Try FFT on each things
	 * 	Do FFT on Taw_polar tensor
	 * Taw_polar <- FFT(Taw_polar)
	 */

	polar_field();	// It performs aforementioned step 1 and step 2.


	/*
	 * step 3	Calculate Taw_polar <--- (-Green_op : Taw_polar)
	 */
	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t k=0;k<N[2];k++){
				double xi2_v[] = {kx_v[i],ky_v[j],kz_v[k]};
				double xi2_v_2 = xi2_v[0]*xi2_v[0]+xi2_v[1]*xi2_v[1]+xi2_v[2]*xi2_v[2];
				if (xi2_v_2>0.0e0){// Only do it for where there is nonzero wave vector
					for (size_t ii=0;ii<dim;ii++)
						for (size_t jj=0;jj<dim;jj++){
							F_hat[in5.g(i,j,k,ii,jj)] = 0.0;
							for (size_t kk=0;kk<dim;kk++)
								for (size_t ll=0;ll<dim;ll++){
									F_hat[in5.g(i,j,k,ii,jj)] += -Taw_polar[in5.g(i,j,k,kk,ll)]*static_cast<complex<double> >(inv_accu_h[in5.g(i,j,k,ii,kk)]*xi2_v[jj]*xi2_v[ll]);
								}
						}
					/*
					 * Actually, it is wrong code
					complex<double> ctmp[dim][dim];
					for (size_t a=0;a<dim;a++)
						for (size_t b=0;b<dim;b++){
							ctmp[a][b]=0.0;
							for (size_t c=0;c<dim;c++)
								for (size_t d=0;d<dim;d++){
									ctmp[a][b] += -static_cast<complex<double> >(green_h[in7.g(i,j,kk,a,b,c,d)])*Taw_polar[in5.g(i,j,kk,c,d)];
								}
						}
					for (size_t a=0;a<dim;a++)
						for (size_t b=0;b<dim;b++){
							Taw_polar[in5.g(i,j,kk,a,b)] = ctmp[a][b];
						}
					 */
				}
				else{	//	Where there is zero wave number, F_avg should be assigned.
					for (size_t a=0;a<dim;a++)
						for (size_t b=0;b<dim;b++){
							//							cout << "Fstep is given as"<<endl;
							//							helpers::print_matrix(F_step,dim,dim);
							F_hat[in5.g(i,j,k,a,b)] = static_cast<complex<double> >(F_step[a*dim+b]*nxnynz);
						}
				}
			}
	/*	step 4 calculate back and same F.
	 * 4.1	Load Taw_polar[k][l] values onto compl_tmp
	 * 4.2	Do ifft
	 * 4.3	Save values on F[k][l]
	 */

	for (size_t k=0;k<dim;k++)
		for (size_t l=0;l<dim;l++){
			//	4.1
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++)
					for (size_t kk=0;kk<N[2];kk++){
						compl_tmp[in3.g(i,j,kk)]=F_hat[in5.g(i,j,kk,k,l)];
					}
			//	4.2
			fft::fftw(N,dim,&compl_tmp[0],&compl_tmp[0],0,nprocs);
			//	4.3
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++)
					for (size_t kk=0;kk<N[2];kk++){
						//	Take real value of the output
						//	My fftw automatically normalizes the output.
						F[in5.g(i,j,kk,k,l)]=(compl_tmp[in3.g(i,j,kk)].real());
					}
		}
}

/*	Convergence is measured based on F and F_prev
 * Average of Frobeneous norm is measured
 *
 */
unsigned int FI_simul::check_conv(double tol){

	size_t N1[]={N[0],N[1],N[2],dim,dim};
	helpers::indexing in5(N1,5);	// Only pass first 5
	size_t nxnynz= N[0]*N[1]*N[2];
	double error = 0.0;

	for (size_t i=0;i<nxnynz;i++){
		double frobenius = 0.0;
		for (size_t k=0;k<dim;k++)
			for (size_t l=0;l<dim;l++){
				frobenius += pow(F[i]-F_prev[i],2);
			}
		error += sqrt(frobenius);
	}
	error /= nxnynz;
	if (error <= tol){
		cout << "FI scheme converged "<< error<<endl;
		//	Essential things that should be done after convergence is attained.
		calc_disp();
		if (!enable_nlocal) save_damage();
		convert_P_2_sig();
		calc_P_F_avg();
		return 1;
	}
	else{
		cout << "FI scheme did not converged "<< error<<endl;
		return 0;
	}


}


/*
 * find_ref_props: Find reference material's property
 * It is defined as maximum parameters of elasticity over a RVE
 */
void FI_simul::find_ref_pros(){


	vector<double> props_list(ela_prop.size(),0);
	//	i: material index, j: property index
	//	Try to find maximum j-property among ela_prop.size() number of materials.
	for (size_t j=0;j<num_ref_props;j++){
		for (size_t i=0;i<num_materials;i++){
			props_list[i]=ela_prop[i][j];
			//			cout << props_list[i]<<endl;
		}
		ref_props[j]=*std::max_element(props_list.begin(),props_list.end())*1.2;
	}
	//	double mu0=ref_props[0]; double lam0 = ref_props[1];

}

/*	After convergence is attained, displacement field will be calculated
 *
 */
void FI_simul::calc_disp(){

	size_t N1[]={N[0],N[1],N[2],dim,dim};
	helpers::indexing in3(N1,3);	// Only pass first 3
	helpers::indexing in4(N1,4);	// Only pass first 4
	helpers::indexing in5(N1,5);	// Only pass first 5

	/*	Step 1: Calculate Perturbation (Polarization) field again
	 *
	 */

	//	Step 1.2. Take FFT of Taw_polar

	polar_field();


	/*	Step 2: calculate disp_h(m) = Acoustic tensor(mk)*i_hat * T_polar(k,L)*kvect[L]
	 *	It is implemented to calculate displacement DOF one by one
	 *	complex_tmp = acoustic tensor(m,k)*i_hat*T_polar(k,L)*kvect(L)
	 *	dummy indices are k and L.
	 */
	complex<double>	i_hat(0.0,1.0);
	for (size_t m=0;m<dim;m++){
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					//	Initialize...
					compl_tmp[in3.g(i,j,kk)] = {0.0,0.0};
					double kvect_v[] = {kx_v[i],ky_v[j],kz_v[kk]};
					//	Only do the calculation for non-zero wave vector grid point.
					double abs_k = kvect_v[0]*kvect_v[0]+kvect_v[1]*kvect_v[1]+kvect_v[2]*kvect_v[2];
					if (abs_k > 0.0){
						for (size_t k=0;k<dim;k++)
							for (size_t L=0;L<dim;L++){
								compl_tmp[in3.g(i,j,kk)] += i_hat*Taw_polar[in5.g(i,j,kk,k,L)]*static_cast<complex<double> >(inv_accu_h[in5.g(i,j,kk,m,k)]*kvect_v[L]);
								//compl_tmp[in3.g(i,j,kk)] += i_hat*Taw_polar[in5.g(i,j,kk,k,L)]*static_cast<complex<double> >(inv_accu_h[in5.g(i,j,kk,k,m)]*kvect_v[L]);
							}
					}
					else{
						//	Set average value of displacement * N[0]N[1][N[2]. Just do nothing
					}
				}
		/*	Step 3: Do Inverse FFT and retrieve displacement field.
		 *
		 */
		fft::fftw(N,dim,&compl_tmp[0],&compl_tmp[0],0,1);
		for (size_t i=0;i<N[0];i++)
			for (size_t j=0;j<N[1];j++)
				for (size_t kk=0;kk<N[2];kk++){
					//	Displacement array is always 3D.
					size_t disp_indx=i*N[1]*N[2]*3+j*N[2]*3+kk*3+m;
					disp[disp_indx] = compl_tmp[in3.g(i,j,kk)].real();
				}
	}

}

/*	init_damage
 *  Build damage model vector
 */
void FI_simul::initalize_damage(unsigned int* mat_indexes){

	size_t N1[]={N[0],N[1],N[2],dim,dim};
	helpers::indexing in3(N1,3);	// Only pass first 3
	helpers::indexing in5(N1,5);	// Only pass first 5

	double F33[3][3];
	math_op::initiali_arr(&F33[0][0],0.0e0,9);
	F33[2][2] = 1.0;
	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t kk=0;kk<N[2];kk++){
				unsigned int indx = mat_indexes[in3.g(i,j,kk)];

				//				double arg_elaprops[]= {ela_prop[indx][0], ela_prop[indx][1]};
				//				double arg_damprops[]= {dam_prop[indx][0], dam_prop[indx][1], dam_prop[indx][2]};
				if (dim==2){
					for (size_t k=0;k<dim;k++)
						for (size_t l=0;l<dim;l++)	{
							F33[k][l]= F[in5.g(i,j,kk,k,l)];
						}
					if (enable_nlocal){
						n_dam.emplace_back('N',2,0,ela_prop[indx],dam_prop[indx],0.0,0.0,0.0,&F33[0][0]);
					}
					else{
						grid_dam.emplace_back('N',2,0,ela_prop[indx],dam_prop[indx],0.0,0.0,&F33[0][0]);
					}
				}
				else{	// dim == 3
					if (enable_nlocal){
						n_dam.emplace_back('N',2,0,ela_prop[indx],dam_prop[indx],0.0,0.0,0.0,&F[in5.g(i,j,kk,0,0)]);
					}
					else{
						grid_dam.emplace_back('N',2,0,ela_prop[indx],dam_prop[indx],0.0,0.0,&F[in5.g(i,j,kk,0,0)]);
					}
				}
			}
	if (enable_nlocal){
		if(debug)	cout << "start initializing vector variables for solving nonlocal damage problem."<<endl;
		size_t nxnynz = N[0]*N[1]*N[2];
		dA12.resize(nxnynz,0.0);
		dA21.resize(nxnynz,0.0);
		dA22.resize(nxnynz,0.0);
		r_phid.resize(nxnynz,0.0);
//		r_nloc.resize(nxnynz,0.0);
		r_nloc=VectorXd::Zero(nxnynz);
		eigfdop.reserve(1);
		eigfdop.emplace_back(nxnynz,nxnynz);
		savefdop.resize(nxnynz,0.0);	//	Only save diagonal term
		FD::build_CD_matrix(eigfdop[0],N,dx,dam_prop[0][4],0);
		if (debug) cout << "building CD matrix is done"<<endl;
		for (size_t i=0;i<nxnynz;i++) savefdop[i] = eigfdop[0].coeffRef(i,i);
		if(debug)cout<<"show save_fdop"<<endl;
		if(debug)helpers::print_matrix(&savefdop[0],1,nxnynz);
		//		fdop('R',1,nxnynz,nxnynz5);
	}

}


void FI_simul::summary(){
	cout<<endl; for(size_t i=0;i<40;i++) cout << "*";	cout<<endl;
	cout <<"Summary of Fixed-Point iteration information"<<endl;
	cout << "Dimension of the problem [" << dim<<"]" << endl;
	cout << "Grid points " << N[0]<<" "<<N[1]<<" "<<N[2]<< endl;
	cout << "Number of processors " << nprocs << endl;
	cout << "preconditioning " << (preconditioning ? "applied" : "not used")  << endl;
	if (preconditioning==1) cout<<"Diagonal preconditioning"<<endl;
	else if (preconditioning==2) cout<<"Incomplete LU preconditioning"<<endl;
	cout << ((restartflag==1) ? "Restart method ":"Not restart method") << endl;
	cout << "number of materials: ["<<num_materials<<"]"<<endl;
	cout << "number of elastic properties (reference material properties): ["<<num_ref_props<<"]"<<endl;
	cout << "number of damage properties (reference material properties): ["<<num_dam_props<<"]"<<endl;
	cout << "Print out reference material's properties "<<"[";
	for(size_t i=0;i<num_ref_props;i++){
		cout<<ref_props[i]; if (i!=num_ref_props-1) cout<<", ";
	}	cout <<"]"<< endl;


	for(size_t i=0;i<40;i++) cout << "*";
	cout<<endl;
}

/*
 * Initialize field quantities using arg_vtk
 * Field quantities to be read are
 * F: deformation gradient
 * P: first piola-Kirchhoff
 * sig: Cauchy stress
 * grid_dam_v: damage
 * grid_ndam_v: kappa
 * displacement: displacement
 */
void FI_simul::restart(string arg_vtk){
	cout << "Restart VTK file"<<endl;
	restartflag=1;
	std::ifstream ifs(arg_vtk);
	string linetxt,dummy,dummy2;
	size_t a,b,c;
	size_t rest_flags[]={0,0,0,0,0,0,0};	//	rest_flags for checking validation of file and field variables
	size_t N2[]={N[0],N[1],N[2],dim,dim};
	helpers::indexing in3(N2,3);	// Only pass first 3
	helpers::indexing in5(N2,5);	// Only pass first 5
	while (std::getline(ifs,linetxt)){
		//		cout << linetxt<<endl;
		if (linetxt.compare(0,10,"DIMENSIONS")==0){
			cout << "Found dimensions..."<<endl;
			std::stringstream iss(linetxt);
			iss >> dummy >> a >> b >> c;
			//	Do verification of grid file.
			if (((N[0]==a)&&(N[1]==b)&&(N[2]==c))==0){
				cout << "Dimension mismatch between Eta file and VTK file for restart"<<endl;
				cout << "N1,N2,N3 "<< N[0]<<" "<<N[1]<<" "<<N[2]<<endl;
				cout << "VTK dimension is..." << a << " "<<b<<" "<<c<<endl;
				exit(-1);
			}
			cout << "dimension is..." << a << " "<<b<<" "<<c<<endl;
			rest_flags[0]=1;
		}
		if (rest_flags[1]==0){	//	Try to find out P
			std::stringstream iss(linetxt);
			iss >> dummy >> dummy2 ;
			if (dummy2.compare("P")==0){//	P found
				cout << "Found P"<<endl;
				getline(ifs,dummy);	//	Skip one line line
				for (size_t k=0;k<N[2];k++)
					for (size_t j=0;j<N[1];j++)
						for (size_t i=0;i<N[0];i++){
							getline(ifs,linetxt);	//	Skip one line
							//							cout << "index i,j,k "<<i<<" "<<j<<" "<<k<<endl;
							std::stringstream remove_space(linetxt);
							for (size_t ii=0;ii<dim;ii++)
								for (size_t jj=0;jj<dim;jj++){
									remove_space >> P[in5.g(i,j,k,ii,jj)];
									//									cout << P[in5.g(i,j,k,ii,jj)]<<" ";
									//									if ((ii==dim-1)&&(jj==dim-1)) cout << endl;
								}
						}
				rest_flags[1]=1;
			}
		}
		if (rest_flags[2]==0){	//	Try to find out Sig
			std::stringstream iss(linetxt);
			iss >> dummy >> dummy2 ;
			if (dummy2.compare("Sig")==0){//	P found
				cout << "Found Sig"<<endl;
				getline(ifs,dummy);	//	Skip one line line
				for (size_t k=0;k<N[2];k++)
					for (size_t j=0;j<N[1];j++)
						for (size_t i=0;i<N[0];i++){
							getline(ifs,linetxt);	//	Skip one line
							//							cout << "index i,j,k "<<i<<" "<<j<<" "<<k<<endl;
							std::stringstream remove_space(linetxt);
							for (size_t ii=0;ii<dim;ii++)
								for (size_t jj=0;jj<dim;jj++){
									remove_space >> sig[in5.g(i,j,k,ii,jj)];
									//									cout << sig[in5.g(i,j,k,ii,jj)]<<" ";
									//									if ((ii==dim-1)&&(jj==dim-1)) cout << endl;
								}
						}
				rest_flags[2]=1;
			}
		}
		if (rest_flags[3]==0){	//	Try to find out damage
			std::stringstream iss(linetxt);
			iss >> dummy >> dummy2 ;
			if (dummy2.compare("damage")==0){//	P found
				cout << "Found damage"<<endl;
				getline(ifs,dummy);	//	Skip one line line
				for (size_t k=0;k<N[2];k++)
					for (size_t j=0;j<N[1];j++)
						for (size_t i=0;i<N[0];i++){
							getline(ifs,linetxt);	//	Skip one line
							//							cout << "index i,j,k "<<i<<" "<<j<<" "<<k<<endl;
							std::stringstream remove_space(linetxt);
							size_t index = in3.g(i,j,k);
							remove_space >> grid_dam_v[index];
							if (grid_dam_v[index] > 1.0e-7) flags[0] = 1;	// Damage initiated flag.
							if(enable_nlocal){
								n_dam[index].dam = grid_dam_v[index];
								n_dam[index].xid = grid_dam_v[index];
							}
							else{
								grid_dam[index].dam = grid_dam_v[index];
								grid_dam[index].xid = grid_dam_v[index];
							}
						}
				rest_flags[3]=1;
			}
		}
		if (enable_nlocal){	//	Try to find kapp only if nonlocal damage model is used.
			if (rest_flags[4]==0){	//	Try to find out nonlocal damage
				std::stringstream iss(linetxt);
				iss >> dummy >> dummy2 ;
				if (dummy2.compare("kappa")==0){//	P found
					cout << "Found kappa"<<endl;
					getline(ifs,dummy);	//	Skip one line line
					for (size_t k=0;k<N[2];k++)
						for (size_t j=0;j<N[1];j++)
							for (size_t i=0;i<N[0];i++){
								getline(ifs,linetxt);	//	Skip one line
								std::stringstream remove_space(linetxt);
								size_t index = in3.g(i,j,k);
								remove_space >> grid_ndam_v[index];
								n_dam[index].ndam = grid_ndam_v[index];
								//								cout << "index i,j,k "<<i<<" "<<j<<" "<<k<<" ndam "<<n_dam[index].ndam<<endl;
							}
					rest_flags[4]=1;
					//					exit(-1);
				}
			}
		}
			if (rest_flags[5]==0){	//	Try to find out F
				std::stringstream iss(linetxt);
				iss >> dummy >> dummy2 ;
				if (dummy2.compare("F")==0){//	P found
					cout << "Found F"<<endl;
					getline(ifs,dummy);	//	Skip one line line
					for (size_t k=0;k<N[2];k++)
						for (size_t j=0;j<N[1];j++)
							for (size_t i=0;i<N[0];i++){
								getline(ifs,linetxt);	//	Skip one line
								//							cout << "index i,j,k "<<i<<" "<<j<<" "<<k<<endl;
								std::stringstream remove_space(linetxt);
								for (size_t ii=0;ii<dim;ii++)
									for (size_t jj=0;jj<dim;jj++){
										remove_space >> F[in5.g(i,j,k,ii,jj)];
										//										cout << F[in5.g(i,j,k,ii,jj)]<<" ";
										//										if ((ii==dim-1)&&(jj==dim-1)) cout << endl;
									}
							}
					rest_flags[5]=1;
					//					exit(-1);
				}
			}
			if (rest_flags[6]==0){	//	Try to find out displacement
				size_t N1[]={N[0],N[1],N[2],3};
				helpers::indexing in4(N1,4);	// Only pass first 3
				std::stringstream iss(linetxt);
				iss >> dummy >> dummy2 ;
				if (dummy2.compare("displacement")==0){//	P found
					cout << "Found displacement"<<endl;
					for (size_t k=0;k<N[2];k++)
						for (size_t j=0;j<N[1];j++)
							for (size_t i=0;i<N[0];i++){
								getline(ifs,linetxt);	//	Skip one line
								//							cout << "index i,j,k "<<i<<" "<<j<<" "<<k<<endl;
								std::stringstream remove_space(linetxt);
								for (size_t ii=0;ii<3;ii++){
									remove_space >> disp[in4.g(i,j,k,ii)];
									//										cout << F[in5.g(i,j,k,ii,jj)]<<" ";
									//										if ((ii==dim-1)&&(jj==dim-1)) cout << endl;
								}
							}
					rest_flags[6]=1;
				}
			}
		}
	ifs.close();
}



FI_simul::~FI_simul() {
	// TODO Auto-generated destructor stub
}

