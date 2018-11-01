#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 16:48:53 2018

@author: jungjaeyong
"""

import os
import sys
file_dir = os.path.dirname(os.path.abspath("mbe_tests.py"))
sys.path.append(file_dir+"/../src") # add relative path
import numpy as np
Dtyp = np.float64

import mesh
import mesh_gen
import MBE

gmesh = mesh.mesh("global")
N=[2,2,2]
L=[1.0,1.0,1.0]
#L=[2.0,2.0,2.0]
#L=[4.0,4.0,4.0]
mesh_gen.structured_hex(N,L,gmesh)

    # MBE tests
tmbe = MBE.MBE(gmesh,mat_type="dense")
# apply boundary conditions
# fix condition
zero_node = gmesh.find_node_indx([0.0,0.0,0.0])
tmbe.bc.add_dirichlet_bc(zero_node,'x',0.0)
tmbe.bc.add_dirichlet_bc(zero_node,'y',0.0)
tmbe.bc.add_dirichlet_bc(zero_node,'z',0.0)
# fix x=0 surface
xnodes=gmesh.find_nodes_on_surf('x',0.0)
tmbe.bc.add_dirichlet_bc(xnodes,'x',0.0)
# stretch x=1 surface
xnodes=gmesh.find_nodes_on_surf('x',1.0)
tmbe.bc.add_dirichlet_bc(xnodes,'x',0.5)
#    tmbe.bc.add_perodic_bc(1,0,[0.0,0.0,0.0])
tmbe.set_activ_dof()
tmbe.assign_d_bc()
for i in range(8):
    tmbe.el_params['hex8'][0][i][0:3]=np.array([0.0,5000.0,7500.0])
tmbe.loop_elem()
tmbe.solve()
tmbe.update_disp()




