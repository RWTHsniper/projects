/*
 * damage.h
 *
 *  Created on: Nov 7, 2017
 *      Author: jaeyong
 */

#ifndef NLOC_DAMAGE_H_
#define NLOC_DAMAGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <vector>
#include "helpers.h"
#include "continuum_mech.h"

class nloc_damage {
private:
	char consti_type;	//	Usually, it is Neo-Hook. 'N'
	unsigned int damfct;	// damage function. Usually, it is 2 and fD = [1-D]^2
	unsigned int hardening;	// index of hardening function. 0: exponential hardening.
	vector<double> ela_prop;	//	Mu, Lam
	vector<double> dam_prop;	//	[Y0,r,s,H,alph]
	double	fD[3];	//	Damage function and its derivatives
	double	F[3][3]; 	// Deformation gradient.
	double	C[3][3]; 	// right-cauchy green tensor.
	double	detC;
	double	trC;
	double	lndetC;
	double	TCF[2];	//	Thermodynamic conjugate forces. Y and qd;

//	Private methods
	void calc_fD();
	void calc_TCF();
	void calc_phid();

public:
	double	phid;
	double	H;
	double	alph;
	double	kloc;
	double	dam;	//	damage state.
	double	ndam;	//	Nonlocal damage term.
	double	xid;	//	damage hardening state.

//	public methods
	nloc_damage(char arg_type,unsigned int arg_damfunct,unsigned int arg_hardfunct,vector<double> &arg_elaprops,vector<double> &arg_damprops,double arg_prevD,double arg_prevND,double arg_prevxid,double* arg_F);	//	Nonlocal damage model
	virtual ~nloc_damage();
	void updateF(double* argF);
	int update_dam();
	double	get_dam(){return dam;};
	double	get_ndam(){return ndam;};
	double	get_xid(){return xid;};
	double	get_fD(){return	fD[0];};
	double	get_phid(){return phid;};
};

#endif /* NLOC_DAMAGE_H_ */
