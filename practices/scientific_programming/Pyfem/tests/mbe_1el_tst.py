#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 11:06:32 2018

@author: localadmin
"""

import os
import sys
file_dir = os.path.dirname(os.path.abspath("mbe_tests.py"))
sys.path.append(file_dir+"/../src") # add relative path
import numpy as np
Dtyp = np.float64
Ityp = np.int64
import mesh
import mesh_gen

import MBE
import postp

"""
Notes on E and poisson's ratio

E= mu*(3*lam+2*mu)/(lam+mu)
pratio = lam/(2*(lam+mu))

mu=7500
lam=5000
E=18000.0
pr = 0.2
"""

# test gauss-quadrature
import GQ
gq_tst = GQ.GQ_hex8()
gq_tst.check()


gmesh = mesh.mesh("global")
N=[2,2,2]
#    N=[5,2,2]
#L=[1.0,1.0,1.0]
L=[2.0,2.0,2.0]
#L=[4.0,4.0,4.0]
mesh_gen.hex8_ref(gmesh)
#gmesh.nodes[3],gmesh.nodes[2] = (gmesh.nodes[2],gmesh.nodes[3])
#gmesh.nodes[6],gmesh.nodes[7] = gmesh.nodes[7],gmesh.nodes[6]
#gmesh.hex8 = np.array([0,1,2,3,4,5,6,7],dtype=Ityp)
#    elem_coords = np.zeros((8,3),dtype=Dtyp)
#    gmesh.elem_node_coords(gmesh.hex8[0],elem_coords)
#    
#    tst_coords= np.array([[-1,1,1],[-1,-1,1],[-1,-1,-1],[-1,1,-1],[1,1,1],[1,-1,1],[1,-1,-1],[1,1,-1]],dtype=Dtyp)
#    
#    gq_tst.calc_J(tst_coords)
#    for i in range(8):
#        print(gq_tst.J[i])
#        
#    gq_tst.calc_dNdXs(tst_coords)

mu=7500
lam=5000
E=18000.0
pr = 0.2

# MBE tests
tmbe = MBE.MBE(gmesh,mat_type="dense",nlgeom=True)
#tmbe = MBE.MBE(gmesh,mat_type="dense",nlgeom=False)
# apply boundary conditions
# fix conditions
# fix x=0 surface in x direction
xnodes=gmesh.find_nodes_on_surf('x',-1.0)
tmbe.bc.add_dirichlet_bc(xnodes,'x',0.0)
# fix y=0 surface in y direction
xnodes=gmesh.find_nodes_on_surf('y',-1.0)
tmbe.bc.add_dirichlet_bc(xnodes,'y',0.0)
# fix z=0 surface in y direction
xnodes=gmesh.find_nodes_on_surf('z',-1.0)
tmbe.bc.add_dirichlet_bc(xnodes,'z',0.0)

# stretch x=1 surface
xnodes=gmesh.find_nodes_on_surf('x',1.0)
tmbe.bc.add_dirichlet_bc(xnodes,'x',0.2)
#    tmbe.bc.add_perodic_bc(1,0,[0.0,0.0,0.0])
tmbe.set_activ_dof()
#    tmbe.assign_d_bc()
for i in range(tmbe.nhex8):        
    for j in range(8):
        tmbe.el_params['hex8'][i][j][0:3]=np.array([1,lam,mu]) # mat_indx,mu,lam
tmbe.NR_iteration(max_it=80)
mypostp = postp.postp(gmesh)
mypostp.add_point_data(tmbe.disp,'displacement_point')
mypostp.write("1el_disp_tst.vtk")
