/*
 * damage.cpp
 *
 *  Created on: Nov 7, 2017
 *      Author: jaeyong
 */

#define DEBUG 0

#include "nloc_damage.h"

using namespace std;

/*
 * nloc_damage
 * Nonlocal damage model is solved.
 */
nloc_damage::nloc_damage(char arg_type,unsigned int arg_damfunct,unsigned int arg_hardfunct,vector<double> &arg_elaprops,vector<double> &arg_damprops,double arg_prevD,double arg_prevND,double arg_prevxid,double* arg_F) {
	// TODO Auto-generated constructor stub

	consti_type = arg_type;
	damfct = arg_damfunct;
	hardening = arg_hardfunct;
	ela_prop=arg_elaprops;	//	[Mu, Lam]
	dam_prop=arg_damprops;	//	[Y0,r,s]
	dam=arg_prevD;
	ndam=arg_prevND;
	xid=arg_prevxid;
	phid = 0;
	kloc = 0;
	H = dam_prop[3];
	alph = dam_prop[4];
	if (consti_type != 'N'){
		cout << "currently, only 'N' type is supported" << endl;
	}
	if (hardening==0){
	}
	else
	{
		cout << "hardening model other than exponential model is not supported. Aborting.."<<endl;
		abort();
	}

	memcpy(&F[0][0],arg_F,9*sizeof(double));
	finite_strain::right_cauchy_green_T(&C[0][0],&F[0][0]);
	detC = math_op::determinant(&C[0][0],3);
	trC = C[0][0]+C[1][1]+C[2][2];
	lndetC = log(detC);
}


void nloc_damage::updateF(double* argF){
	memcpy(&F[0][0],argF,9*sizeof(double));
	finite_strain::right_cauchy_green_T(&C[0][0],&F[0][0]);
	detC = math_op::determinant(&C[0][0],3);
	trC = C[0][0]+C[1][1]+C[2][2];
	lndetC = log(detC);
}


/*	calc_fD
 * calculate damage function
 */
void nloc_damage::calc_fD(){
	if (damfct == 2){
		fD[0] = (1-dam)*(1-dam);
		fD[1] = 2*(dam-1);
		fD[2] = 2;
	}
	else{
		cout << "Other type of damage function is not supported."<<endl;
		exit(-1);
	}
}

/*	calc_TCF
 * calculate thermodynamic conjugate forces
 * TCF[0]	=	Y;	TCF[1]	=	qd
 */
void nloc_damage::calc_TCF(){
	if (consti_type == 'N'){
		//	  Y=-fD(2)*(mu/2.0d0*(trC-3.0d0-log_detC)+lambda/4.0d0*(detC-1.0d0-log_detC))
		double mu = ela_prop[0];
		double lambda = ela_prop[1];
		double H = dam_prop[3];
		TCF[0] = -fD[1]*(mu/2*(trC-3-lndetC)+lambda/4*(detC-1-lndetC)) - H* (dam-ndam);
	}
	if (hardening==0){
//		  q_d = r * (1.0d0 - dexp(-s*xi_d))
		double r = dam_prop[1];
		double s = dam_prop[2];
		  TCF[1] = r * (1 - exp(-s*xid));
	}
}

/*
 * calculate damage loading function.
 */
void nloc_damage::calc_phid(){
	double Y0 = dam_prop[0];
	calc_fD();
	calc_TCF();
	phid = TCF[0]-(Y0+TCF[1]);
}

/*	update_dam
 * try to update damage
 * If it is one shot clean to pass, number 0 is returned.
 * Otherwise, damage is updated and 1 will be returned.
 */
int nloc_damage::update_dam(){
/*
 *
 * 	    log_detC = log(detC)
dY_dd_lam = -fD(3)*(mu/2.0d0*(trC-3.0d0-log_detC)+lam/4.0d0*(detC-1.0d0-log_detC))-H
dq_dd_lam = d_param(1)*d_param(3)*dexp(-d_param(3)*(xi_d)) ! xi_d: value for the next step
K_loc= dY_dd_lam - dq_dd_lam
 *
 */

	if (DEBUG){
		cout << "update_damage"<<endl;
		cout << "printF"<<endl;
		helpers::print_matrix(&F[0][0],3,3);
		cout << "printC"<<endl;
		helpers::print_matrix(&C[0][0],3,3);
		cout << "print detC, lndetC and trC "<< detC<<" "<<lndetC<<" "<<trC<<endl;
	}

	if (hardening==0){
    	calc_phid();
		if (phid > 1.0e-7){
			double mu = ela_prop[0];
			double lambda = ela_prop[1];
			double r = dam_prop[1];
			double s = dam_prop[2];
			double H = dam_prop[3];
			double dY_ddlam = -fD[2]*(mu/2*(trC-3-lndetC)+lambda/4*(detC-1-lndetC))-H;
		    double dqd_ddlam = r*s*exp(-s*xid);
		    kloc = dY_ddlam - dqd_ddlam;
			return 1;
		}
		else{
//			phid = 0.0;
			return 0;
		}
	}	//	hardening == 0

cout << "Hardening model is invalid."<<endl;
return -1;
}

nloc_damage::~nloc_damage() {
	// TODO Auto-generated destructor stub
}

