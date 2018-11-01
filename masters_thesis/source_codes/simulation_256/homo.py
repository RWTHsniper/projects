#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:48:36 2017

@author: jaeyong
"""

import os
import re
import numpy as np
import operator
import vtk
from vtk.util import numpy_support as VN
#from os.path import isfile, join


class file_reader:
    def __init__(self, title, dim):
        if (dim == 2):
            data = [[],[]]
        elif (dim == 3):
            data = [[],[],[]]
        lendata = 0
        with open(title) as f:
            for line in f:
                dum_data = line.split()
                lendata = lendata + 1
                for i in range(dim):
                    data[i].append(float(dum_data[i]))
        self.outpdata = np.zeros((lendata,dim),dtype=np.float,order='C')
        for i in range(lendata):
            for j in range(dim):
                self.outpdata[i][j] = data[j][i]	
    def return_data(self):
        return self.outpdata
# Destroyer
	def __del__(self):
		class_name = self.__class__.__name__
		print class_name, "destroyed"

"""
Main section of my code
"""

outpfile = "homo_stre_str.txt"
tol = 1.0e-6
arr = os.listdir('./vtkfiles')
arr.sort()
# Import data
titleoffile = "z_stress_strain.txt"
reader = file_reader(titleoffile,2)
str_strain_data = reader.return_data()
# Find the peack point
maxindx = str_strain_data[:,1].argmax()

"""
Example of VTK reader. Read only one VTK file firstly.
"""
"""
#	VTK Reader initialization
filename = "simul00010.vtk"
reader = vtk.vtkRectilinearGridReader()
reader.SetFileName(filename)
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()

data = reader.GetOutput()
#data.GetDimensions()
#mat_indx = VN.vtk_to_numpy(data.GetPointData().GetArray('mat_indx'))
#P=VN.vtk_to_numpy(data.GetPointData().GetArray('P'))
sig=VN.vtk_to_numpy(data.GetPointData().GetArray('Sig'))
dam=VN.vtk_to_numpy(data.GetPointData().GetArray('damage'))
"""
file = open(outpfile,"w") 
for vtkf in arr:
    vtknumber=map(int,re.findall('\d+',vtkf))
#    print vtknumber
#	Try homogenization after reaching peack
    if (vtknumber>=maxindx):
        filename = vtkf
        reader = vtk.vtkRectilinearGridReader()
        reader.SetFileName("./vtkfiles/"+filename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        dam_flag=VN.vtk_to_numpy(data.GetPointData().GetArray('dam_flag'))
        sig=VN.vtk_to_numpy(data.GetPointData().GetArray('Sig'))
        dam=VN.vtk_to_numpy(data.GetPointData().GetArray('damage'))
#        F=VN.vtk_to_numpy(data.GetPointData().GetArray('F'))
        dims = data.GetDimensions()
        avg_sig = np.array([0.0,0.0,0.0,0.0])
        dimension=dims[0]*dims[1]*dims[2]
        counter = 0
        for i in range(dimension):
            if (dam[i]>tol and dam_flag[i]==1):
                avg_sig += sig[i,:]
                counter += 1
#            print "damage bigger than 0"
#            print "summation of sigma "
        file.writelines(str(float(str_strain_data[vtknumber,0]))+" "+str(avg_sig[0]/counter)+"\n") 
file.close()           









