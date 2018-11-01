#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 11:32:02 2018

@author: localadmin

Gaussian Quadrature info 

nn: number of nodes
ngq: number of gauss points
N: shape functions for each gauss point (self.ngq,self.nn)
N has shape function values at gauss points
w: weight for each gauss point
gp: coordinates of gauss point

"""

import sys
import math
import numpy as np

deftype = np.float64
Dtyp = deftype

class GQ_hex8:
    def __init__(self):
        self.nn = 8
        self.nen = 8
        self.ngq = 8
        self.ndim = 3
        self.sdim = 3
        self.N = np.zeros((self.ngq,self.nn),dtype=deftype)
        self.w = np.ones(self.ngq,dtype=deftype)
        self.gp = np.zeros((self.ngq,self.ndim),dtype=deftype)
        self.dNdxi = np.zeros((self.ngq,self.nn),dtype=deftype)
        self.dNdeta = np.zeros((self.ngq,self.nn),dtype=deftype)
        self.dNdzeta = np.zeros((self.ngq,self.nn),dtype=deftype)
        self.J = np.zeros((self.ngq,self.ndim,self.ndim),dtype=deftype)
        self.detJ = np.zeros((self.ngq),dtype=deftype)
        self.dNdXs = np.zeros((self.ngq,self.ndim,self.nn),dtype=deftype)

        # assign gps
        isqrt3 = 1.0/math.sqrt(3.0)
        self.gp[0] = np.array([-isqrt3,isqrt3,isqrt3],dtype=deftype)
        self.gp[1] = np.array([-isqrt3,-isqrt3,isqrt3],dtype=deftype)
        self.gp[2] = np.array([-isqrt3,-isqrt3,-isqrt3],dtype=deftype)
        self.gp[3] = np.array([-isqrt3,isqrt3,-isqrt3],dtype=deftype)
        self.gp[4] = np.array([isqrt3,isqrt3,isqrt3],dtype=deftype)
        self.gp[5] = np.array([isqrt3,-isqrt3,isqrt3],dtype=deftype)
        self.gp[6] = np.array([isqrt3,-isqrt3,-isqrt3],dtype=deftype)
        self.gp[7] = np.array([isqrt3,isqrt3,-isqrt3],dtype=deftype)

        fct = 8.0

        # calculate shape functions
        for i in range(self.ngq): # i: index of gauss point
            self.N[i][0] = (1.0-self.gp[i][0])*(1.0+self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.N[i][1] = (1.0-self.gp[i][0])*(1.0-self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.N[i][2] = (1.0-self.gp[i][0])*(1.0-self.gp[i][1])*(1.0-self.gp[i][2])/fct
            self.N[i][3] = (1.0-self.gp[i][0])*(1.0+self.gp[i][1])*(1.0-self.gp[i][2])/fct
            self.N[i][4] = (1.0+self.gp[i][0])*(1.0+self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.N[i][5] = (1.0+self.gp[i][0])*(1.0-self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.N[i][6] = (1.0+self.gp[i][0])*(1.0-self.gp[i][1])*(1.0-self.gp[i][2])/fct
            self.N[i][7] = (1.0+self.gp[i][0])*(1.0+self.gp[i][1])*(1.0-self.gp[i][2])/fct

        # calculate derivatives of shape function in the reference element
        for i in range(self.ngq):
            # calculate dNdxi
            self.dNdxi[i][0] = -(1.0+self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.dNdxi[i][1] = -(1.0-self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.dNdxi[i][2] = -(1.0-self.gp[i][1])*(1.0-self.gp[i][2])/fct
            self.dNdxi[i][3] = -(1.0+self.gp[i][1])*(1.0-self.gp[i][2])/fct
            self.dNdxi[i][4] = (1.0+self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.dNdxi[i][5] = (1.0-self.gp[i][1])*(1.0+self.gp[i][2])/fct
            self.dNdxi[i][6] = (1.0-self.gp[i][1])*(1.0-self.gp[i][2])/fct
            self.dNdxi[i][7] = (1.0+self.gp[i][1])*(1.0-self.gp[i][2])/fct
            # calculate dNdeta
            self.dNdeta[i][0] = (1.0-self.gp[i][0])*(1.0+self.gp[i][2])/fct
            self.dNdeta[i][1] = -(1.0-self.gp[i][0])*(1.0+self.gp[i][2])/fct
            self.dNdeta[i][2] = -(1.0-self.gp[i][0])*(1.0-self.gp[i][2])/fct
            self.dNdeta[i][3] = (1.0-self.gp[i][0])*(1.0-self.gp[i][2])/fct
            self.dNdeta[i][4] = (1.0+self.gp[i][0])*(1.0+self.gp[i][2])/fct
            self.dNdeta[i][5] = -(1.0+self.gp[i][0])*(1.0+self.gp[i][2])/fct
            self.dNdeta[i][6] = -(1.0+self.gp[i][0])*(1.0-self.gp[i][2])/fct
            self.dNdeta[i][7] = (1.0+self.gp[i][0])*(1.0-self.gp[i][2])/fct
            # calculate dNdzeta
            self.dNdzeta[i][0] = (1.0-self.gp[i][0])*(1.0+self.gp[i][1])/fct
            self.dNdzeta[i][1] = (1.0-self.gp[i][0])*(1.0-self.gp[i][1])/fct
            self.dNdzeta[i][2] = -(1.0-self.gp[i][0])*(1.0-self.gp[i][1])/fct
            self.dNdzeta[i][3] = -(1.0-self.gp[i][0])*(1.0+self.gp[i][1])/fct
            self.dNdzeta[i][4] = (1.0+self.gp[i][0])*(1.0+self.gp[i][1])/fct
            self.dNdzeta[i][5] = (1.0+self.gp[i][0])*(1.0-self.gp[i][1])/fct
            self.dNdzeta[i][6] = -(1.0+self.gp[i][0])*(1.0-self.gp[i][1])/fct
            self.dNdzeta[i][7] = -(1.0+self.gp[i][0])*(1.0+self.gp[i][1])/fct

# For a given set of coordinates, calculate jacobian
    def calc_J(self,coords):
        for igp in range(self.ngq):
            self.J[igp][0][0]=np.dot(self.dNdxi[igp],coords[:,0])
            self.J[igp][0][1]=np.dot(self.dNdxi[igp],coords[:,1])
            self.J[igp][0][2]=np.dot(self.dNdxi[igp],coords[:,2])
            self.J[igp][1][0]=np.dot(self.dNdeta[igp],coords[:,0])
            self.J[igp][1][1]=np.dot(self.dNdeta[igp],coords[:,1])
            self.J[igp][1][2]=np.dot(self.dNdeta[igp],coords[:,2])
            self.J[igp][2][0]=np.dot(self.dNdzeta[igp],coords[:,0])
            self.J[igp][2][1]=np.dot(self.dNdzeta[igp],coords[:,1])
            self.J[igp][2][2]=np.dot(self.dNdzeta[igp],coords[:,2])
            self.detJ[igp] = np.linalg.det(self.J[igp])
#   The same thing can be done in the below manner
#        j1 = np.array([self.dNdxi[igp],self.dNdeta[igp],self.dNdzeta[igp]],dtype=deftype)
#        np.matmul(j1,coords,J)

# calculate dNdx, dNdy, dNdz
# [dNdX]          [dNdxi]
# [dNdY] = invJ * [dNdeta]
# [dNdZ]          [dNdzeta] 
    def calc_dNdXs(self,coords):
        self.calc_J(coords)
        for igq in range(self.ngq):
            invJ = np.linalg.inv(self.J[igq])
            np.matmul(invJ,np.array([self.dNdxi[igq],self.dNdeta[igq],self.dNdzeta[igq]],dtype=Dtyp),self.dNdXs[igq])

# use natural coordinate to find its N.
    def calc_N(self,nat_coord):
        fct = 8.0
        Nv = np.zeros((self.nn),dtype=deftype)
        Nv[0] = (1.0-nat_coord[0])*(1.0+nat_coord[1])*(1.0+nat_coord[2])/fct
        Nv[1] = (1.0-nat_coord[0])*(1.0-nat_coord[1])*(1.0+nat_coord[2])/fct
        Nv[2] = (1.0-nat_coord[0])*(1.0-nat_coord[1])*(1.0-nat_coord[2])/fct
        Nv[3] = (1.0-nat_coord[0])*(1.0+nat_coord[1])*(1.0-nat_coord[2])/fct
        Nv[4] = (1.0+nat_coord[0])*(1.0+nat_coord[1])*(1.0+nat_coord[2])/fct
        Nv[5] = (1.0+nat_coord[0])*(1.0-nat_coord[1])*(1.0+nat_coord[2])/fct
        Nv[6] = (1.0+nat_coord[0])*(1.0-nat_coord[1])*(1.0-nat_coord[2])/fct
        Nv[7] = (1.0+nat_coord[0])*(1.0+nat_coord[1])*(1.0-nat_coord[2])/fct
        return Nv        

#   This is used to check whether N, dNdxi,dNdeta,dNdzeta are correct
    def check(self):
        tol = 1.0e-6
        delt = 1.0e-7
        fini_dNdxi = np.zeros((self.nn),dtype=deftype)
        print ("Doing shape function summation check")
        for i in range(self.ngq):
            if (math.fabs(1.0-self.N[i].sum())> tol):
                print (str(i)+"-th GP has wrong N = "+str(self.N[i].sum()))
                sys.exit()
        print ("Shape function test passed")
        
        print ("Let's try derivatives of shape functions")
        print ("dNdxi goes on")
        for i in range(self.ngq):
            nat_coord = self.gp[i]
            nat_coord[0] = nat_coord[0] + delt
            Nvr = self.calc_N(nat_coord)
            nat_coord[0] = nat_coord[0] - 2.0*delt
            Nvl = self.calc_N(nat_coord)
            fini_dNdxi = (Nvr-Nvl)/(2.0*delt)
            diff = self.dNdxi[i]-fini_dNdxi
            if (np.linalg.norm(diff) > tol):
                print("Derivatives of shape functions are incorrect at "+str(i)+"-th GP")
                print(fini_dNdxi)
                print(self.dNdxi[i])
                sys.exit()

        print ("dNdeta goes on")
        for i in range(self.ngq):
            nat_coord = self.gp[i]
            nat_coord[1] = nat_coord[1] + delt
            Nvr = self.calc_N(nat_coord)
            nat_coord[1] = nat_coord[1] - 2.0*delt
            Nvl = self.calc_N(nat_coord)
            fini_dNdxi = (Nvr-Nvl)/(2.0*delt)
            diff = self.dNdeta[i]-fini_dNdxi
            if (np.linalg.norm(diff) > tol):
                print("Derivatives of shape functions are incorrect at "+str(i)+"-th GP")
                print(fini_dNdxi)
                print(self.dNdeta[i])
                sys.exit()

        print ("dNdzeta goes on")
        for i in range(self.ngq):
            nat_coord = self.gp[i]
            nat_coord[2] = nat_coord[2] + delt
            Nvr = self.calc_N(nat_coord)
            nat_coord[2] = nat_coord[2] - 2.0*delt
            Nvl = self.calc_N(nat_coord)
            fini_dNdxi = (Nvr-Nvl)/(2.0*delt)
            diff = self.dNdzeta[i]-fini_dNdxi
            if (np.linalg.norm(diff) > tol):
                print("Derivatives of shape functions are incorrect at "+str(i)+"-th GP")
                print(fini_dNdxi)
                print(self.dNdzeta[i])
                sys.exit()
        print("Derivative tests passed!")
                        

    def __del__(self):
        class_name = self.__class__.__name__
        print (class_name, "destroyed")










