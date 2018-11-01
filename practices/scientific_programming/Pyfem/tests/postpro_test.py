#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 14:59:50 2018

@author: jungjaeyong
"""



import os
import sys
file_dir = os.path.dirname(os.path.abspath("postpro_test.py"))
sys.path.append(file_dir+"/../src") # add relative path
import numpy as np
from mesh import *
from mesh_gen import *
from postp import *
from prepro import *

Ityp = np.int64

"""
def main(argv):
    # My code here
    N=[2,2,2] 
    L=[1.0,1.0,1.0]
    hexa =structured_hex(0,N,L)

if __name__ == "__main__":
    main(sys.argv)
    
    """
    # My code here
N=[3,3,3]
#N=[8,8,8]
L=[1.0,1.0,1.0]
gmesh = mesh("global")
structured_hex(N,L,gmesh)
mypostp = postp(gmesh)

## test postprocessing
# test cell data
tensor_celldata = np.zeros([gmesh.nel,3,3],dtype=np.float64)
for i in range(tensor_celldata.shape[0]):
    for j in range(3):
        for k in range(3):
            tensor_celldata[i][j][k]=i+j*3.0+k

scalar_celldata = np.zeros([gmesh.nel,1],dtype=np.float64)
for i in range(gmesh.nel):
    scalar_celldata[i] = i

vector_celldata = np.zeros([gmesh.nel,3],dtype=np.float64)
for i in range(gmesh.nel):
    vector_celldata[i,:] = np.array([i,i,i],dtype=np.float64)    
    
            
mypostp.add_cell_data(tensor_celldata,'tensor_test')
for i in range(tensor_celldata.shape[0]):
    for j in range(3):
        for k in range(3):
            tensor_celldata[i][j][k] = 2

# create point data
vector_point = np.zeros([gmesh.nn,3],dtype=np.float64)
factor = 0.01
for i in range(gmesh.nn):
    vector_point[i] = np.array([factor*gmesh.nodes[i][0],factor*gmesh.nodes[i][1],factor*gmesh.nodes[i][2]],dtype=np.float64)

mypostp.add_cell_data(tensor_celldata,'tensor_test')
mypostp.add_cell_data(vector_celldata,'vector_cell')
mypostp.add_cell_data(scalar_celldata,'scalar_cell_indx')
mypostp.add_point_data(vector_point,'vector_point')
mypostp.write("test.vtk")


