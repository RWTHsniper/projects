#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 21:56:53 2018

@author: jungjaeyong
"""

import vtk
import numpy as np
from vtk.util import numpy_support


class postp:
    def __init__(self,mesh):
        # Add the points and hexahedron to an unstructured grid
        self.uGrid = vtk.vtkUnstructuredGrid()
        # make object of points
        self.points = vtk.vtkPoints()
#        self.points.SetNumberOfPoints(mesh.nodes.shape[0])
        for i in range(mesh.nodes.shape[0]):
            self.points.InsertNextPoint([mesh.nodes[i,0],mesh.nodes[i,1],mesh.nodes[i,2]])
        self.uGrid.SetPoints(self.points)
            
#        if (hasattr(mesh,"hex8")):
#        indx_h8 = [mesh.elements.index(x) for x in mesh.elements if x[0]=='hex8']         # index of 'hex8'
        if (mesh.elements['hex8'][1] > 0):
            # Add the hexahedron to a cell array
            for i in range(mesh.hex8.shape[0]):
                hexa8 = vtk.vtkHexahedron()
                for j in range(8):
                    hexa8.GetPointIds().SetId(j,mesh.hex8[i,j])
                self.uGrid.InsertNextCell(hexa8.GetCellType(),hexa8.GetPointIds())
    
    def add_cell_data(self,data,dat_name):
        if (data.dtype==np.float64):
            VTKtype = vtk.VTK_DOUBLE
        if (data.shape[1]==1):    # scalar or vector quantity
            VTKarray = numpy_support.numpy_to_vtk(num_array=data,deep=True,array_type=VTKtype)
            self.uGrid.GetCellData().SetScalars(VTKarray)
        elif ((data.shape[1]==3) and (len(data.shape)==2)):     # scalar or vector quantity
            VTKarray = numpy_support.numpy_to_vtk(num_array=data,deep=True,array_type=VTKtype)
            self.uGrid.GetCellData().SetVectors(VTKarray)            
        elif (len(data.shape)==3):  # tensor quantity
            VTKarray = numpy_support.numpy_to_vtk(num_array=data.reshape(data.shape[0],data.shape[1]*data.shape[2]),deep=True,array_type=VTKtype)
            self.uGrid.GetCellData().SetTensors(VTKarray)
        VTKarray.SetName(dat_name)

    def add_point_data(self,data,dat_name):
        if (data.dtype==np.float64):
            VTKtype = vtk.VTK_DOUBLE
        if (data.shape[1]==1):    # scalar or vector quantity
            VTKarray = numpy_support.numpy_to_vtk(num_array=data,deep=True,array_type=VTKtype)
            self.uGrid.GetPointData().SetScalars(VTKarray)
        elif ((data.shape[1]==3) and (len(data.shape)==2)):     # scalar or vector quantity
            VTKarray = numpy_support.numpy_to_vtk(num_array=data,deep=True,array_type=VTKtype)
            self.uGrid.GetPointData().SetVectors(VTKarray)            
        elif (len(data.shape)==3):  # tensor quantity
            VTKarray = numpy_support.numpy_to_vtk(num_array=data.reshape(data.shape[0],data.shape[1]*data.shape[2]),deep=True,array_type=VTKtype)
            self.uGrid.GetPointData().SetTensors(VTKarray)
        VTKarray.SetName(dat_name)
        
    def write(self,filename):
        #Write the grid in a vtu file
        modelwriter = vtk.vtkUnstructuredGridWriter()
        modelwriter.SetFileName(filename)
        modelwriter.SetInputData(self.uGrid)
        modelwriter.Write()
        
    def __del__(self):
        class_name = self.__class__.__name__
        print (class_name, "destroyed")
        

    
"""    
    def put_tensor_cell_data(self,data,dat_name):
        if (data.dtype==np.float64):
            VTKtype = vtk.VTK_DOUBLE
        VTKarray = numpy_support.numpy_to_vtk(num_array=data.reshape(data.shape[0],data.shape[1]*data.shape[2]),deep=True,array_type=VTKtype)
        VTKarray.SetName(dat_name)
        self.uGrid.GetCellData().SetScalars(VTKarray)
#        self.uGrid.SetActiveScalars

    def put_scalar_cell_data(self,data,dat_name):
        if (data.dtype==np.float64):
            VTKtype = vtk.VTK_DOUBLE
        VTKarray = numpy_support.numpy_to_vtk(num_array=data,deep=True,array_type=VTKtype)
        VTKarray.SetName(dat_name)
        self.uGrid.GetCellData().SetTensors(VTKarray)
"""

            

