
/*
 ============================================================================
 Name        : res_Eig_FFT_2D_vidyasagar.c
 Author      :
 Version     :
 Copyright   : Your copyright notice
 Description : Hello OpenMP World in C
 ============================================================================
 */

#define debug 0

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cstdlib>	//	used to run system() command.
#include <ctime>	//	timer
#include <ctype.h>	//	isdigit

//	Classes and functions written by me
//#include "ffts.h"
#include "continuum_mech.h"
#include "grid2d.h"
#include "vtks.h"
#include "FIsimul.h"
#include "damage.h"
#include "nloc_damage.h"

using namespace std;

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

//	restart, restartfile should be specified
	size_t restart = 1;	//	Restart option
	string restartfile = "./vtkfiles/simul00772.vtk";	//	Restart VTK file
	std::ofstream outp;
	size_t int_extrac = 0;
	size_t initial_loading_step=0;
	if (restart){//	restart on/off
		std::ifstream inp("z_stress_strain.txt");
		outp.open("r_z_stress_strain.txt");
		string extracted_num = "";
		for (size_t i=0;i<restartfile.size();i++){
			if (isdigit(restartfile[i])){
				extracted_num+=restartfile[i];
			}
		}
		cout << "extracted_num = "<<extracted_num<<endl;
		std::stringstream iss(extracted_num);
		 iss >> int_extrac;
		cout << "extracted_num = "<<int_extrac<<endl;
		for (size_t i=0;i<int_extrac+1;i++){
			string linetxt;
			std::getline(inp,linetxt);
			outp << linetxt<<endl;
		}
		inp.close();
		system("rm -r r_vtkfiles");
		const int dir_err = system("mkdir -p r_vtkfiles");
		if (-1 == dir_err)
		{
			printf("Error creating r_vtkfiles directory!n");
			exit(1);
		}
		initial_loading_step = int_extrac+1;	// Make it next loading step by adding 1.
	}
	else{
		system("rm -r vtkfiles");
		const int dir_err = system("mkdir -p vtkfiles");
		if (-1 == dir_err)
		{
			printf("Error creating vtkfiles directory!n");
			exit(1);
		}
		outp.open("z_stress_strain.txt");
	}


	/*	Control paramter boards.
	 *
	 */

	/*	Step 1
	 *
	 */

	size_t N[]={256,256,1};	//	Define Nx,Ny,Nz
	string gridtitle="CEta1_256";	char gridtype='C';

//	size_t N[]={128,128,1};	//	Define Nx,Ny,Nz
//	size_t N[]={64,64,1};	//	Define Nx,Ny,Nz
//	size_t N[]={32,32,1};	//	Define Nx,Ny,Nz
//	size_t N[]={16,16,1};	//	Define Nx,Ny,Nz
//	size_t N[]={8,8,1};	//	Define Nx,Ny,Nz
//	size_t N[]={4,4,1};	//	Define Nx,Ny,Nz
//	size_t N[]={2,2,1};	//	Define Nx,Ny,Nz
	double	L[]	={100e-3,100e-3,100e-3};
	/*	step 2	assign grid
	 *
	 */
//	string gridtitle="CEta1_128";	char gridtype='C';
//	string gridtitle="CEta1_64";	char gridtype='C';
//	string gridtitle="CEta1_32";	char gridtype='C';
//	string gridtitle="CEta1_16";	char gridtype='C';
//	string gridtitle="CEta1_8";	char gridtype='C';
//	string gridtitle="CEta1_4";	char gridtype='C';
	unsigned int makehomog	=	0;
	unsigned int nprocs = 1;
	/*	step 3
	 *	initialize fixed point scheme
	 */
	unsigned int nlocdam = 1;
	unsigned int preconditioner = 1;	//	Preconditioner for CG scheme.
	unsigned int ela_dam_num[] = {2,5};	// number of elastic and damage parameters
	unsigned int control_flags[]	={nprocs,nlocdam,preconditioner};
//	double alph = 5.0e-7, H = 1.0e4;
//	double alph = 1.0e-7, H = 1.0e4;
//	double alph = 5.0e-8, H = 1.0e4;
	double alph = 2.5e-8, H = 1.0e4;
//	double alph = 1.0e-8, H = 1.0e4;
	//	Don't forget to initilize FI.mat_list, FI.ela_prop and FI.dam_prop
	/*	step 4
	 * loadstep define
	 */
	unsigned int n_load = 3500;
	unsigned int max_iter = 2500;	//	maximum number of iteration for Fixed point scheme
//	double incr = 1.0e-4;
	double incr = 5.0e-5;
//	double incr = 1.0e-5;
	double tol = 1.0e-7;	//	0.5 % tolerance w.r.t. increment


//	unsigned int drawarr[] = {10,20};
	unsigned int drawarr[] = {5,10};
//	unsigned int drawarr[] = {2,4};
	unsigned int draw = drawarr[0];


	/*	simulation steps
*/
/*	Step 1	Define grid space to define domain of the simulation
	 *
	 */
	grid_2d grid(N[0],N[1],N[2],L[0],L[1],L[2]);
	/*	step 2	Assign grid
	 * By default, the grid is initialized as homogeneous one with index 0.
	 */
	if (!makehomog){
		grid.read_grid(gridtitle,gridtype);
	}

	grid.summary();

	/*	Step 3	Initialize Fixed-Point scheme
	 * Load useful info to the FIsimul class including grid information and material properties
	 */
	//	Firstly, I will run elastic RVE.
	FI_simul FI(grid,grid.num_materials,ela_dam_num[0],ela_dam_num[1],control_flags);
	size_t sdim = FI.dim;
	if(debug)	cout<<"Initialization of FI_simul is done."	<<endl;
	//	mat_list and ela_prop should be defined.
	FI.mat_list[0]=0;
	FI.ela_prop[0][0]=	7500.0; FI.ela_prop[0][1]=5000.0;	//	[Mu,Lambda]
	FI.dam_prop[0][0] = 5.0; FI.dam_prop[0][1] = 50.0; FI.dam_prop[0][2] = 0.5;	//	[Y0,r,s]
	FI.dam_prop[0][3] = H; FI.dam_prop[0][4] = alph;//	[H,alph]

	if (!makehomog){
		FI.mat_list[1]=1;
		double factorM1 = 1.05;
		FI.ela_prop[1][0]=	factorM1*7500.0; FI.ela_prop[1][1]=factorM1*5000.0;	//	[Mu,Lambda]
		FI.dam_prop[1][0] = 1.0e12; FI.dam_prop[1][1] = 50.0; FI.dam_prop[1][2] = 0.5;	// [Y0,r,s]
		FI.dam_prop[1][3] = H; FI.dam_prop[1][4] = alph;//	[H,alph]
	}

	if(debug) cout << "damage module is going to be initialized"<<endl;
	FI.initalize_damage(grid.get_mat_indx());	//	damage model should be inisitlized in FI_simul object.
	if(debug) cout << "damage module initialized"<<endl;
	FI.find_ref_pros();	// Find out reference material's properties
	FI.def_accu_green();	//	Define Inverse of accustic tensor and accustic tensor
	//	Initialize Stress field
	FI.calc_P_field(grid);
	if(restart) FI.restart(restartfile);
	FI.summary();	// Shows up summary of FI scheme

	/*
	 * Check Whether VTK file is correctly read
	 * After validation, this section can be killed.
	 */
	/*
	//	Pad leading zeros to make it read in order in Paraview
	string leadingzero_i=helpers::put_leading_zeros(5,int_extrac);
	string vtktitle;
	if (restart) vtktitle = "./r_vtkfiles/simul"+leadingzero_i;
	else vtktitle = "./vtkfiles/simul"+leadingzero_i;
	//	Initialize vtk object
	vtks vtk(grid.N,grid.dx);
	vtk.init_r_grid_vtkfile(vtktitle);	// initialize vtk file headers
	vtk.init_r_grid_point_data();	// start writing point data in vtk file.
	vtk.write_scalar_point_data(grid.get_mat_indx(),static_cast<std::string>("mat_indx"),1);
	vtk.write_scalar_point_data(FI.get_P(),static_cast<std::string>("P"),sdim*sdim);
	vtk.write_scalar_point_data(FI.get_F(),static_cast<std::string>("F"),sdim*sdim);
	vtk.write_scalar_point_data(FI.get_sig(),static_cast<std::string>("Sig"),sdim*sdim);
	vtk.write_vector_point_data(FI.get_disp(),"displacement");	// It is writing vector...
	vtk.write_scalar_point_data(FI.get_damage(),static_cast<std::string>("damage"),1);
	if (FI.enable_nlocal) vtk.write_scalar_point_data(FI.get_nloc_damage(),static_cast<std::string>("kappa"),1);
	vtk.close();
*/

	//	Set initial F
	double F_step[sdim][sdim];
	if (restart){	//	Initial F of restart.
		FI.calc_P_F_avg();
		for (size_t i=0;i<sdim;i++){
			for (size_t j=0;j<sdim;j++){
				F_step[i][j]=FI.Favg[i*sdim+j];
			}
		}
//		F_step[0][0]=1.00055;	If you want to make the starting step very exact, you should put this line but it is not very necessary.
	}
	else{
		for (size_t i=0;i<sdim;i++){
			for (size_t j=0;j<sdim;j++){
				F_step[i][j]=0.0e0;F_step[j][i]=0.0e0;
			}
		}
		F_step[0][0] = 1.0e0;
	}

	/*
	 * Step 4 Start loading steps
	 * 4.1	prepare parameters for loading
	 */
	//	4.1
	//	P-stretch, cauchy stress curve.

	//	Start Load steps
	 clock_t startTime = clock(); //Start timer
	  time_t rawtime;	//	local time
	  time (&rawtime);
	  printf ("The current local time is: %s", ctime (&rawtime));


	for	(size_t i=initial_loading_step;i<n_load;i++){	// int_extrac: start load steop
		//	Define Favg on your own.
		if (FI.flags[0]){	//	If it has damage evolution, try smaller timestep.
			F_step[0][0] += incr/(static_cast<double>(drawarr[1])/static_cast<double>(drawarr[0]));
			draw = drawarr[1];
		}
		else{
			F_step[0][0] += incr;	//	Purely elastic state, try bigger step size.
		}
		F_step[1][1] = 1.0e0/(F_step[0][0]); F_step[0][1]=0.0;F_step[1][0]=0.0;
		if (sdim==3){
			F_step[1][1] = 1.0e0/sqrt(F_step[0][0]);
			F_step[2][2] = F_step[1][1];
		}
		/*	Print loadstep info
		 *
		 */
		for (size_t i=0;i<20;i++) cout<<"*";
		cout << endl<<"runtime "<<	((clock()-startTime)/CLOCKS_PER_SEC)/60 <<" mins "<< ((clock()-startTime)/CLOCKS_PER_SEC)%60 << " secs";
		cout <<endl<<"Loadstep info, "<<i<<endl;
		cout <<endl<<"F_grad"<<endl;
		helpers::print_matrix(&F_step[0][0],FI.dim,FI.dim);

		for (size_t j=0;j<max_iter;j++){	// Number of fixed point iteration
			FI.run_simul(&F_step[0][0]);
			FI.umat(grid);

			unsigned int conv = FI.check_conv(tol);
			if (conv == 0){
				if (j==max_iter-1){
					cout << "Fixed point iteration scheme reached maximum number of iteration. aborting..."<<endl;
					abort();
				}
				continue;
			}
			else{
				cout << "Solution converged!!"<<endl;

				//	Prepare saving data on z_stress_strain.txt file.
				double sig[3][3], F33[3][3],P33[3][3],C33[3][3];
				math_op::initiali_arr(&F33[0][0],0.0e0,9);
				math_op::initiali_arr(&P33[0][0],0.0e0,9);
				for (size_t i=0;i<sdim;i++)
					for (size_t j=0;j<sdim;j++)	{
						F33[i][j]=F_step[i][j];
						P33[i][j]=FI.Pavg[i*sdim+j];
					}
				F33[2][2]=1.0e0;P33[2][2]=0.0e0;
				finite_strain::right_cauchy_green_T(&C33[0][0],&F33[0][0]);
				elasticity::convert_P_2_Sig(&sig[0][0],&P33[0][0],&F33[0][0]);
				double	lam_max =math_op::eig_power(&C33[0][0],3);
				lam_max = sqrt(lam_max);
				outp <<  lam_max<<" "<<sig[0][0]<<endl;
				//				outp <<  F_step[0][0]<<" "<<sig[0][0]<<endl;

				//	Write down vtk file
				if ((i % draw)==0){
					//	Pad leading zeros to make it read in order in Paraview
					string leadingzero_i=helpers::put_leading_zeros(5,i);
					string vtktitle;
					if (restart) vtktitle = "./r_vtkfiles/simul"+leadingzero_i;
					else vtktitle = "./vtkfiles/simul"+leadingzero_i;
					//	Initialize vtk object
					vtks vtk(grid.N,grid.dx);
					vtk.init_r_grid_vtkfile(vtktitle);	// initialize vtk file headers
					vtk.init_r_grid_point_data();	// start writing point data in vtk file.
					vtk.write_scalar_point_data(grid.get_mat_indx(),static_cast<std::string>("mat_indx"),1);
					vtk.write_scalar_point_data(FI.get_dam_flag(),static_cast<std::string>("dam_flag"),1);
					vtk.write_scalar_point_data(FI.get_P(),static_cast<std::string>("P"),sdim*sdim);
					vtk.write_scalar_point_data(FI.get_F(),static_cast<std::string>("F"),sdim*sdim);
					vtk.write_scalar_point_data(FI.get_sig(),static_cast<std::string>("Sig"),sdim*sdim);
					vtk.write_vector_point_data(FI.get_disp(),"displacement");	// It is writing vector...
					vtk.write_scalar_point_data(FI.get_damage(),static_cast<std::string>("damage"),1);
					if (FI.enable_nlocal) vtk.write_scalar_point_data(FI.get_nloc_damage(),static_cast<std::string>("kappa"),1);
					vtk.close();
				}
/*
				cout << "Print F33"<< endl;
				helpers::print_matrix(&F33[0][0],3,3);
				cout << "Print P.avg"<<endl;
				helpers::print_matrix(&FI.Pavg[0],2,2);
				cout << "Print P33"<<endl;
				helpers::print_matrix(&P33[0][0],3,3);

				cout << "Print sig avg"<< endl;
				helpers::print_matrix(&sig[0][0],3,3);
*/
				break;




				/*	I will start from working on this...
				 *
				 */
			}
		}
	}
	outp.close();

	/*	Test elasticity module firstly
	 *
	 */
	/*
	double l11 = 2.0e0;
	double Ftest[9]={l11,0.0,0.0,0.0,1.0/l11,0.0,0.0,0.0,1.0};
	double Ptest[9],Sigtest[9], Stest[9];
	cout << "Test module*******************"<<endl;
	cout<< "Print F"<<endl; helpers::print_matrix(Ftest,3,3);
	elasticity::hyper::neo_hook::piola_kirch_stress(Ptest,Ftest,7500.0,5000.0,1);
	cout<< "Print P"<<endl; helpers::print_matrix(Ptest,3,3);
	elasticity::convert_P_2_Sig(&Sigtest[0],&Ptest[0],&Ftest[0]);
	elasticity::hyper::neo_hook::piola_kirch_stress(Stest,Ftest,7500.0,5000.0,2);
	cout<< "Print S"<<endl; helpers::print_matrix(Stest,3,3);
	cout<< "Print Sig"<<endl; helpers::print_matrix(Sigtest,3,3);
	 */

	/*
	 * Jobs that I should do.
	 * 1.	Make P field visualized using paraview
	 * 2.	Then, make stress-strain curve.
	 * 3.	Implement elastic-damage model.
	 */


	/*	complex type test
	 *
	 */


	cout << "End program"<<endl;
	return 0;
}


