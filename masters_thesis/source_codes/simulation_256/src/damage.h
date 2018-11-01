/*
 * damage.h
 *
 *  Created on: Nov 7, 2017
 *      Author: jaeyong
 */

#ifndef DAMAGE_H_
#define DAMAGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <vector>
#include "helpers.h"
#include "continuum_mech.h"

class damage {
private:
	char consti_type;	//	Usually, it is Neo-Hook. 'N'
	size_t damfct;	// damage function. Usually, it is 2 and fD = [1-D]^2
	size_t hardening;	// index of hardening function. 0: exponential hardening.
	vector<double> ela_prop;	//	Mu, Lam
	vector<double> dam_prop;	//	[Y0,r,s]
	double	fD[3];	//	Damage function and its derivatives
	double	F[3][3]; 	// Deformation gradient.
	double	C[3][3]; 	// right-cauchy green tensor.
	double	detC;
	double	trC;
	double	lndetC;
	double	H;	//	Penalty term for nonlocal damage
	double	alph;
	double	TCF[2];	//	Thermodynamic conjugate forces. Y and qd;

//	Private methods
	void calc_fD();
	void calc_TCF();
	void calc_phid();

public:
	double	dam;	//	Initial damage state.
	double	xid;	//	Initial damage hardening state.
	double	phid;

//	public methods
	damage(char arg_type,size_t arg_damfunct,size_t arg_hardfunct,vector<double> &arg_elaprops,vector<double> &arg_damprops,double arg_prevD,double arg_prevxid,double* arg_F);
	virtual ~damage();
	void updateF(double* argF);
	int update_dam();
	double	get_dam(){return dam;};
	double	get_xid(){return xid;};
	double	get_fD(){return	fD[0];};
};

#endif /* DAMAGE_H_ */
