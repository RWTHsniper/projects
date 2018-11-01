#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 23:23:25 2018

@author: jungjaeyong

Material models for MBE

"""

"""
inp
   C: C=F'*F to be calculated
   F: deformation gradient
"""        

import numpy as np
Dtyp = np.float64
import math

def r_cauchy_green_T(F,C=None):
    if (C is None):
        return np.matmul(F.transpose(),F)
    else:
        np.matmul(F.transpose(),F,C[:,:])   
           

def green_lagrangian_T(C,E=None):
    if (E is None):
        return 0.5*(C-np.identity(C.shape[0]))
    else:
        E[:,:] = 0.5*(C-np.identity(C.shape[0]))

def T2_Mandel(T2):
    sdim = T2.shape[0]
    M = np.zeros((sdim*2,1),dtype=Dtyp)
    sqrt2 = math.sqrt(2)
    for i in range(sdim):
        M[i] = T2[i,i]
        M[i+sdim] = sqrt2*T2[(i)%sdim,(i+1)%sdim]
    return M

def Mandel2_T(M,T=None):
    isqrt2 = 1.0/math.sqrt(2)
    if len(M)==6:
        T = np.zeros((3,3),dtype=Dtyp)
        T[0,0]=M[0]
        T[1,1]=M[1]
        T[2,2]=M[2]
        T[0,1]=M[3]*isqrt2
        T[1,0]=T[0,1]
        T[1,2]=M[4]*isqrt2
        T[2,1]=T[1,2]
        T[2,0]=M[5]*isqrt2
        T[0,2]=T[2,0]
        return T
    if (T is not None):
        if len(M)==6:
            T[0,0]=M[0]
            T[1,1]=M[1]
            T[2,2]=M[2]
            T[0,1]=M[3]*isqrt2
            T[1,0]=T[0,1]
            T[1,2]=M[4]*isqrt2
            T[2,1]=T[1,2]
            T[2,0]=M[5]*isqrt2
            T[0,2]=T[2,0]
        

"""
small strain tensor calculator
"""
def small_strain(F,eps=None):
    if (eps is None):
        eps = 0.5*(F+F.transpose()) - np.identity(F.shape[0])
        return eps
    else:
        eps[:,:] = 0.5*(F+F.transpose()) - np.identity(F.shape[0])

"""
lin_elasticity
dsde: dsig/deps
F: deformation gradient
S: 2nd Piola-Kirchhoff stress


"""


def lin_elasticity(dsde,F,S,params,**kwargs):   
    lam = params[0]
    mu = params[1]
    eps = small_strain(F)
    sdim = F.shape[0]
    for i in range(sdim):
        dsde[i,i] = lam + 2.0*mu
        dsde[i+sdim,i+sdim] = 2.0*mu
        for j in range(i+1,sdim):
            dsde[i,j] = lam
            dsde[j,i] = lam
    Meps = T2_Mandel(eps) # eps in Mandel notation
    sig = Mandel2_T(np.matmul(dsde,Meps))
    invF = np.linalg.inv(F)
    S[:,:] = np.linalg.det(F) * np.matmul(np.matmul(invF,sig),invF.transpose())



"""
D: algorithmic tangent
F: deformation gradient
S: 2nd Piola-Kirchhoff stress
params: [mu,lam]
"""


def Neo_Hook(dSdE,F,S,params,**kwargs):
    Dtyp = S.dtype
    
    def calc_stress(mu,lam,C,S):
        #   Calculate second PK first
        invC = np.linalg.inv(C)
        S[:,:] = mu * (np.identity(S.shape[0])-invC)+lam/2.0*(np.linalg.det(C)-1.0)*invC

#        invC = np.linalg.inv(C)
#        detC = np.linalg.det(C)
#        for i in range(3):
#            for j in range(3):
#                # firstly, start with calculating the second PK stress
#                S[i,j] = mu*(1.0*(i==j)-invC[i,j])+lam/2.0*(detC-1.0)*invC[i,j]       
                
    lam = params[0]
    mu = params[1]
    C = np.zeros((3,3),dtype=Dtyp)
    Cp = np.zeros((3,3),dtype=Dtyp)
    r_cauchy_green_T(F,C)
    calc_stress(mu,lam,C,S)
    St = np.zeros((3,3),dtype=Dtyp)
    Stmp = np.zeros((3,3),dtype=Dtyp)    
    calc_stress(mu,lam,C,St)
    if (kwargs is not None):
        for key, value in kwargs.items():
            if (key == "finite_diff"):
                if (value == True):
                    # continue writing finite difference method
#                    pertb = min(C[C.nonzero()])*1.0e-3
                    factor = 1.0e-7
                    minp = 1.0e-10
                    # dSdC 00
                    Cp[:,:] = C[:,:]
                    pertb = np.max([minp,factor * Cp[0,0]])
                    Cp[0,0] = Cp[0,0]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[0,0] = (Stmp[0,0]-St[0,0])/pertb
                    Cp[0,0] = C[0,0]
                    # dSdC 11
                    pertb = np.max([minp,factor * Cp[1,1]])
                    Cp[1,1] = Cp[1,1]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[1,1] = (Stmp[1,1]-St[1,1])/pertb
                    Cp[1,1] = C[1,1]
                    # dSdC 22
                    pertb = np.max([minp,factor * Cp[2,2]])
                    Cp[2,2] = Cp[2,2]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[2,2] = (Stmp[2,2]-St[2,2])/pertb
                    Cp[2,2] = C[2,2]
                    # dSdC 01
                    pertb = np.max([minp,factor * Cp[1,1]])
                    Cp[1,1] = Cp[1,1]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[0,1] = (Stmp[0,0]-St[0,0])/pertb
                    Cp[1,1] = C[1,1]
                    # dSdC 12
                    pertb = np.max([minp,factor * Cp[2,2]])
                    Cp[2,2] = Cp[2,2]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[1,2] = (Stmp[1,1]-St[1,1])/pertb
                    Cp[2,2] = C[2,2]
                    # dSdC 20
                    pertb = np.max([minp,factor * Cp[0,0]])
                    Cp[0,0] = Cp[0,0]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[2,0] = (Stmp[2,2]-St[2,2])/pertb
                    Cp[0,0] = C[0,0]
                    # symmetry
                    dSdE[1,0]=dSdE[0,1]
                    dSdE[2,1]=dSdE[1,2]
                    dSdE[0,2]=dSdE[2,0]
                    # shear... 
                    # dSdC 33 (0,1)(0,1)
                    pertb = np.max([minp,factor * Cp[0,1]])                    
                    Cp[0,1] = Cp[0,1]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[3,3] = (Stmp[0,1]-St[0,1])/(2.0*pertb)
                    Cp[0,1] = C[0,1]
                    # dSdC 44 (1,2)(1,2)
                    pertb = np.max([minp,factor * Cp[1,2]])
                    Cp[1,2] = Cp[1,2]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[4,4] = (Stmp[1,2]-St[1,2])/(2.0*pertb)
                    Cp[1,2] = C[1,2]
                    # dSdC 55 (2,0)(2,0)
                    pertb = np.max([minp,factor * Cp[2,0]])
                    Cp[2,0] = Cp[2,0]+pertb
                    calc_stress(mu,lam,Cp,Stmp)
                    dSdE[5,5] = (Stmp[2,0]-St[2,0])/(2.0*pertb)
                    Cp[2,0] = C[2,0]
                    # convert dSdC into dSdE. dSdE = 2*dSdC
                    dSdE[:,:]= 2.0*dSdE

"""
Hooke material model is simply used.
c11 = lam + 2*mu
c12 = lam
c44 = mu (I believe c44 should be 2*mu. Let's use it in a conventional manner.)
e.g.
lam = 60.41e9
mu=23.17e9
c11 = 106.75e9
c12 = 60.41e9
c44(mu) = 23.17e9
"""
def hook(dSdE,F,S,params,**kwargs):
    c11 = params[0]
    c12 = params[1]
    c44 = params[2]
    E = green_lagrangian_T(r_cauchy_green_T(F))
    sdim = F.shape[0]
    for i in range(sdim):
        dSdE[i,i] = c11
        dSdE[i+sdim,i+sdim] = 2.0*c44
        for j in range(i+1,sdim):
            dSdE[i,j] = c12
            dSdE[j,i] = c12
    ME = T2_Mandel(E) # eps in Mandel notation
    S[:,:] = Mandel2_T(np.matmul(dSdE,ME))
                    
    
"""
dC = 1e-6*norm(C);
% Ignore what dF dPdF P stands for
% Calculate dSdC firstly.
F = C;
dF = dC;
dPdF = zeros(6,6);
P = calc_S(mu,lam,C);
% 11
F2 = F + dF*[1 0 0; 0 0 0; 0 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(1,1) = (P2(1,1)-P(1,1))/dF;
% 22
F2 = F + dF*[0 0 0; 0 1 0; 0 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(2,2) = (P2(2,2)-P(2,2))/dF;
% 33
F2 = F + dF*[0 0 0; 0 0 0; 0 0 1];
P2 = calc_S(mu,lam,F2);
dPdF(3,3) = (P2(3,3)-P(3,3))/dF;
% 12
F2 = F + dF*[0 0 0; 0 1 0; 0 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(1,2) = (P2(1,1)-P(1,1))/dF;
% 13
F2 = F + dF*[0 0 0; 0 0 0; 0 0 1];
P2 = calc_S(mu,lam,F2);
dPdF(1,3) = (P2(1,1)-P(1,1))/dF;
% 21
F2 = F + dF*[1 0 0; 0 0 0; 0 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(2,1) = (P2(2,2)-P(2,2))/dF;
% 23
F2 = F + dF*[0 0 0; 0 0 0; 0 0 1];
P2 = calc_S(mu,lam,F2);
dPdF(2,3) = (P2(2,2)-P(2,2))/dF;
% 31
F2 = F + dF*[1 0 0; 0 0 0; 0 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(3,1) = (P2(3,3)-P(3,3))/dF;
% 32
F2 = F + dF*[0 0 0; 0 1 0; 0 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(3,2) = (P2(3,3)-P(3,3))/dF;
% 44 (1,2)/(1,2)
F2 = F + dF*[0 1 0; 1 0 0; 0 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(4,4) = (P2(1,2)-P(1,2))/(2*dF);
% 55 (2,3)/(2,3)
F2 = F + dF*[0 0 0; 0 0 1; 0 1 0];
P2 = calc_S(mu,lam,F2);
dPdF(5,5) = (P2(2,3)-P(2,3))/(2*dF);
% 66 (3,1)/(3,1)
F2 = F + dF*[0 0 1; 0 0 0; 1 0 0];
P2 = calc_S(mu,lam,F2);
dPdF(6,6) = (P2(3,1)-P(3,1))/(2*dF);

dPdF = 2*dPdF;  % Convert dSdC into dSdE

<<<<<<< HEAD
"""                
