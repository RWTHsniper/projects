#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 14:55:01 2018

@author: jungjaeyong


input
N: array of the number of nodes
L: length of an RVE

output
hex_elems

""" 

import numpy as np
Dtyp = np.float64
uItyp = np.uint64

"""
Generate reference element of a hex8
"""
def hex8_ref(mesh):
    mesh.hex8 = np.array([[0,1,2,3,4,5,6,7]],dtype=uItyp)
    mesh.nodes = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],dtype=Dtyp)  
    #update number of nodes and elements
    mesh.update_nel_nn()

class structured_hex:
    def __init__(self,N,L,mesh):
        Nn = N[0]*N[1]*N[2]
        Nx=N[0]; Ny=N[1]; Nz=N[2];
        Ne = (N[0]-1)*(N[1]-1)*(N[2]-1)
        dx = [L[i]/(N[i]-1.0) for i in range(3)]
        # set nodes
        mesh.nodes = np.zeros([Nn,3],dtype=Dtyp)
        cnt = 0
        for k in range(Nz):
            for j in range(Ny):
                for i in range(Nx):
                    mesh.nodes[cnt,:]=[i*dx[0],j*dx[1],k*dx[2]]
                    cnt += 1
                    
        # assign hexa mesh
        # index of 'hex8'
        if (mesh.elements['hex8'][1]==0):    # "Hex8" is not found
            mesh.hex8 = np.zeros([Ne,8],dtype=uItyp)
            counter = 0
            for k in range(Nz-1):
                for j in range(Ny-1):
                    for i in range(Nx-1):
                        ind = i+j*Nx+k*Nx*Ny
                        mesh.hex8[counter,:] = [ind,ind+1,ind+Nx+1,ind+Nx,ind+Nx*Ny,ind+1+Nx*Ny,ind+Nx+1+Nx*Ny,ind+Nx+Nx*Ny]
                        counter+=1
        else:
            for k in range(Nz-1):
                for j in range(Ny-1):
                    for i in range(Nx-1):
                        ind = i+j*Nx+k*Nx*Ny
                        mesh.hex8 = np.append(mesh.hex8,[[ind,ind+1,ind+Nx+1,ind+Nx,ind+Nx*Ny,ind+1+Nx*Ny,ind+Nx+1+Nx*Ny,ind+Nx+Nx*Ny]],axis=0)
                        counter+=1
                        
        #update number of nodes and elements
        mesh.update_nel_nn()
            
    def __del__(self):
        class_name = self.__class__.__name__
        print (class_name, "destroyed")

        
