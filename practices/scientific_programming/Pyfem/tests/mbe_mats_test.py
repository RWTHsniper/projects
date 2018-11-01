#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 13:43:18 2018

@author: localadmin
"""

import MBE_mats
import numpy as np
Dtyp = np.float64

lam=5000.0
mu=7500.0
params = [lam,mu]


F = np.array([[1.01,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],dtype=Dtyp)
S = np.zeros((3,3),dtype=Dtyp)
sig = np.zeros((3,3),dtype=Dtyp)
dsigdeps = np.zeros((6,6),dtype=Dtyp)
dSdE = np.zeros((6,6),dtype=Dtyp)
MBE_mats.lin_elasticity(dsigdeps,F,sig,params)   
MBE_mats.Neo_Hook(dSdE,F,S,params,finite_diff=True)

np.set_printoptions(precision=2)
print("Difference in stress")
print(S-sig)
print("Difference in algorithmic tangent")
print(dSdE-dsigdeps)

