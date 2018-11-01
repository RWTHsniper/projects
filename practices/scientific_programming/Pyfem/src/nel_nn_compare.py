#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 23:48:50 2018

@author: jungjaeyong
"""

import mesh
import mesh_gen
import os
import sys
import numpy as np
Dtyp = np.float64

file_dir = os.path.dirname(os.path.abspath("nel_nn_compare.py"))
sys.path.append(file_dir+"/../src") # add relative path

N=[64,64,64]
#N=[8,8,8]
L=[1.0,1.0,1.0]
gmesh = mesh.mesh("global")
mesh_gen.structured_hex(N,L,gmesh)


print("compare whether it is better to use global matrix or mesh-free method")

elem_matr = np.zeros((gmesh.nel,24,24),dtype=Dtyp)
ndof = gmesh.nn*3
global_matr = np.zeros((ndof,ndof),dtype=Dtyp)

if (sys.getsizeof(elem_matr)> sys.getsizeof(global_matr)):
    print("element matrices consume more memory " + str(sys.getsizeof(elem_matr))+" "+str(sys.getsizeof(global_matr)))
else:
    print("global matrix takes up more memory " + str(sys.getsizeof(elem_matr))+" "+str(sys.getsizeof(global_matr)))

print ("global/elem_matrixs = "+str(sys.getsizeof(global_matr)/sys.getsizeof(elem_matr)))

