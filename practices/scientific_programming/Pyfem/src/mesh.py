#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 08:40:33 2018

@author: jungjaeyong

Elements are...
Hex8
Tet4
elements = {'dictionary': (index of first element, the number of elements)}

"""

import numpy as np
Dtyp = np.float64
Ityp = np.int64
import math

class mesh:
    def __init__(self,name):
        self.name = name
        self.nodes = np.zeros([1,3],dtype=np.float64)
        self.elements = {'hex8':[0,0],'tet4':[0,0]}
        self.nn = 0 # total number of nodes
        self.nel = 0 # total number of elements
        self.hex8 = np.zeros([0,8],dtype=np.float64)
        self.tet4 = np.zeros([0,4],dtype=np.float64)

    def update_nel_nn(self):        
        self.nn = self.nodes.shape[0]
        # assign element information
        self.elements['hex8'][1] = self.hex8.shape[0]
        self.elements['tet4'][1] = self.tet4.shape[0]
        tmp_list = list(self.elements)
        for i in range(len(tmp_list)):
            for j in range(i+1,len(tmp_list)):
                self.elements[tmp_list[j]][0] = self.elements[tmp_list[j]][0] + self.elements[tmp_list[i]][1]

        for inel in self.elements:
            self.nel = self.nel + self.elements[inel][1]   # update total number of elements
            
    def find_node_indx(self,pos):
        tol = 1.0e-6
        save_residual = 1.0
        save_indx = 0
        for indx,node in enumerate(self.nodes):
            diff = node-pos
            err = (diff.dot(diff))**(0.5)
            if err < save_residual:
                save_indx = indx
                save_residual = err
        if (save_residual > tol):
            save_indx = -1            
        return save_indx

#   find nodes on a plane. eg) nodes on x=0.
#   in this case, plane = 'x', val = 0
#   returns list of nodes
    def find_nodes_on_surf(self,plane,val):
        tol = 1.0e-6
        outp_lists = []
        if plane == 'x':
            for indx,node in enumerate(self.nodes):
                diff = math.fabs(node[0]-val)
                if (diff < tol):
                    outp_lists.append(indx)
        if plane == 'y':
            for indx,node in enumerate(self.nodes):
                diff = math.fabs(node[1]-val)
                if (diff < tol):
                    outp_lists.append(indx)
        if plane == 'z':
            for indx,node in enumerate(self.nodes):
                diff = math.fabs(node[2]-val)
                if (diff < tol):
                    outp_lists.append(indx)
        return np.array(outp_lists,dtype=Ityp)
                            
    def elem_node_coords(self,elem_arr,elem_coords):
        for i in range(len(elem_arr)):
            elem_coords[i]=self.nodes[elem_arr[i]]



    def __del__(self):
        class_name = self.__class__.__name__
        print (class_name, "destroyed")

    
        
        
        
        
        
        