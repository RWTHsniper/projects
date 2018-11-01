/*
 * FIsimul.h
 *
 *  Created on: Nov 2, 2017
 *      Author: jaeyong
 */

#ifndef FISIMUL_H_
#define FISIMUL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>     // std::cout
#include <sstream>
#include <string>
#include <algorithm>    // std::min_element, std::max_element
#include <math.h>       /* pow */
#include <fstream>      // std::ofstream
#include "grid2d.h"
#include "vtks.h"
#include "math_op.h"
#include "helpers.h"
#include "continuum_mech.h"
#include "ffts.h"
#include "damage.h"
#include "nloc_damage.h"
#include "indexing.h"
//#include "sparmat.h"
#include <Eigen/Sparse>
#include "fini_diff.h"

//#include "grid2d.h"

using namespace std;
using namespace Eigen;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat; // declares a column-major sparse matrix
class FI_simul {

private:	// private variables
	unsigned int num_materials;
//	grid_2d*	grid;
	size_t N[3];
	double dx[3];
	unsigned int num_ref_props;
	unsigned int num_dam_props;
	vector<double> ref_props;
	unsigned int	nprocs;
	//	double	C0[3][3][3][3];
	vector<double> kx;
	vector<double> ky;
	vector<double> kz;
	vector<double> kx_v;
	vector<double> ky_v;
	vector<double> kz_v;
	size_t restartflag;

	/*	3D case
 * double inv_accu_h[3][3];	//	Inverse of accustic tensor in Frequency domain
	double accu_h[3][3];	//	Accustic tensor in Frequency domain
 *
 */
	vector<size_t> dam_flag;	//	[Nx][Ny][Nz] damage zole flag
	vector<double> inv_accu_h;	//	Inverse of accustic tensor in Frequency domain	[Nx][Ny][2][2]
//	vector<double> accu_h;	//	Accustic tensor in Frequency domain	[Nx][Ny][2][2]
//	vector<double> green_h;	//	[Nx][Ny][2][2][2][2]	To save memory. I will only use inv_accu_h
//	Define variables to be solved
	vector<double> disp;	//	[Nx][Ny][Nz][3]	Always dimension is 3
	vector<double> F;	//	[Nx][Ny][dim]dim]
	vector<double> F_prev;	//	[Nx][Ny][dim]dim]
	vector<double> P;	//	[Nx][Ny][dim][dim]
	vector<double> sig;	//	[Nx][Ny][dim][dim]
	vector<double> grid_dam_v;	//[Nx][Ny][Nz]
//	vector<double> grid_ndam_v;	//[Nx][Ny][Nz]
	VectorXd grid_ndam_v;
	//	vector<spar::spar_mat>	fdop;	//	solve FD-approximated problem.	I made it like this to make it initialized later.
	vector<SpMat> eigfdop; //	solve FD-approximated problem.
	vector<double> savefdop;	//	Stores original value of fdop
	vector<double> dA12;	//	It is obvious that it is always -1 or 0.
	vector<double> dA21;
	vector<double> dA22;
//	vector<double> r_nloc;	//	Residual vector of nonlocal damage field
	VectorXd	r_nloc;
	vector<double> r_phid;	//	Residual vector of phi_d
	ConjugateGradient<SpMat, Lower|Upper,Eigen::IdentityPreconditioner> cgI;
	ConjugateGradient<SpMat, Lower|Upper,Eigen::DiagonalPreconditioner<double>> cgD;
	ConjugateGradient<SpMat, Lower|Upper,Eigen::IncompleteLUT<double>> cgLU;

	vector<std::complex<double> > compl_tmp;
	vector<std::complex<double> > Taw_polar;	//	Polarization field [Nx][Ny][3][3]	We will only use [Nx][Ny][2][2]
	vector<std::complex<double> > F_hat;	//	Polarization field [Nx][Ny][3][3]	We will only use [Nx][Ny][2][2]
	vector<std::complex<double> > coef_2diff;	// coefficient to calculate 2nd order differential term.
	vector<damage> grid_dam;	//[Nx][Ny][Nz]	Local damage
	vector<nloc_damage> n_dam;	//[Nx][Ny][Nz]	Nonlocal damage


//	Private methods
	void def_wave_vects();
	void save_prev_F();
	void umat();
	void polar_field();
	void update_damage();
	void update_nloc_damage(size_t maxiter,double rel_tol);


public:
//	Public variables
	size_t dim;
	//	mat_list, ela_prop, dam_prop is explicitly defined in main function.
	vector<unsigned int>	mat_list;
	vector<vector<double> > ela_prop;	//	[Mu, Lam]
	vector<vector<double> > dam_prop;	//	[Y0,r,s]
	vector<double>	Pavg;
	vector<double>	Favg;
	vector<int>	flags;	//[damage_flag,]
	unsigned int enable_nlocal;	//	enable nonlocal damage model or not;
	unsigned int preconditioning;	//	Whether PCG or CG is used for solving nonlocal damage update.
	/*	For 2D
		double	Pavg[2][2];
		double	Favg[2][2];
	*/

//	Public methods
	FI_simul(grid_2d& arg_grid,unsigned int arg_num_materials,unsigned int num_ela_props,unsigned int num_arg_dam_props,unsigned int* arg_control_flags);
	void find_ref_pros();
	void summary();
	void def_accu_green();
	void calc_P_field(grid_2d& grid);
	void run_simul(double* F_step);
	unsigned int check_conv(double tol);
	void umat(grid_2d& grid);
	void calc_P_F_avg();
	void calc_disp();
	void convert_P_2_sig();
	size_t* get_dam_flag(){return &dam_flag[0];};
	double* get_P(){return &P[0];};
	double* get_F(){return &F[0];};
	double* get_sig(){return &sig[0];};
	double* get_disp(){return &disp[0];};
	double* get_damage(){return &grid_dam_v[0];};
	double* get_nloc_damage(){return &grid_ndam_v[0];};
	void initalize_damage(unsigned int* mat_indexes);
	void save_damage();	//	save damage on continuous memory map for visualization. It is called after convergence is attained.
	void restart(string arg_vtk);

	virtual ~FI_simul();
};

#endif /* FISIMUL_H_ */
