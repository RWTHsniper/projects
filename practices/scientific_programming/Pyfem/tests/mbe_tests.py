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

import mesh
import mesh_gen

import MBE
import postp



def main():

    # test gauss-quadrature
    import GQ
    gq_tst = GQ.GQ_hex8()
    gq_tst.check()


    gmesh = mesh.mesh("global")
#    N=[2,2,2]
#    N=[5,2,2]
#    N=[10,2,2]
#    N=[10,3,3]
    N=[10,10,10]
    L=[1.0,1.0,1.0]
    #L=[2.0,2.0,2.0]
    #L=[4.0,4.0,4.0]
    mesh_gen.structured_hex(N,L,gmesh)
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


    # MBE tests
    tmbe = MBE.MBE(gmesh,mat_type="dense",nlgeom=True)
#    tmbe = MBE.MBE(gmesh,mat_type="dense",nlgeom=False)
    n_step = 20
    stretch_rate = 1.0

    for iterate in range(n_step):
        # apply boundary conditions
        # fix conditions
        # fix x=0 surface in x direction
        xnodes=gmesh.find_nodes_on_surf('x',0.0)
        tmbe.bc.replace_dirichlet_bc(xnodes,'x',0.0)
        # fix y=0 surface in y direction
        xnodes=gmesh.find_nodes_on_surf('y',0.0)
        tmbe.bc.replace_dirichlet_bc(xnodes,'y',0.0)
        # fix z=0 surface in y direction
        xnodes=gmesh.find_nodes_on_surf('z',0.0)
        tmbe.bc.replace_dirichlet_bc(xnodes,'z',0.0)

        # stretch x=1 surface
        xnodes=gmesh.find_nodes_on_surf('x',L[0])
        tmbe.bc.replace_dirichlet_bc(xnodes,'x',stretch_rate*(iterate+1)/n_step)
    #    tmbe.bc.add_perodic_bc(1,0,[0.0,0.0,0.0])
        tmbe.set_activ_dof()
    #    tmbe.assign_d_bc()
        for i in range(tmbe.nhex8):
            for j in range(8):
#                tmbe.el_params['hex8'][i][j][0:3]=np.array([0,5000.0,7500.0]) # Neo-Hooke
                tmbe.el_params['hex8'][i][j][0:4]=np.array([2,106.75e9,60.41e9,23.17e9]) # Hooke in DAMASK
        tmbe.NR_iteration(max_it=30)
        mypostp = postp.postp(gmesh)
        mypostp.add_point_data(tmbe.disp,'displacement_point')
        mypostp.write("disp_tst_"+str(iterate)+".vtk")
#        mypostp.write("disp_NH_"+str(iterate)+".vtk")

if __name__ == "__main__":
    main()

