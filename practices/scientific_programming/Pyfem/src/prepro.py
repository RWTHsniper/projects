#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 20:50:42 2018

@author: jungjaeyong
"""


import vtk
import numpy as np
from vtk.util import numpy_support


class prepro:
    def __init__(self):               
        self.reader = vtk.vtkUnstructuredGridReader()
        self.reader.ReadAllTensorsOn()
        self.reader.ReadAllVectorsOn()
        self.reader.ReadAllScalarsOn()


    def read(self,filename):
        #Write the grid in a vtu file
        self.reader.SetFileName(filename)
        self.reader.Update()     
        self.data = self.reader.GetOutput()

    def assign_geom(self,mesh):
        # initialize mesh.elements
        for eltypes in list(mesh.elements.keys()):
            mesh.elements[eltypes]=[0,0]
        mesh.nn = self.data.GetNumberOfPoints()
        mesh.nodes = np.zeros([mesh.nn,3],dtype=np.float64) # resize mesh.nodes
        for i in range(mesh.nn):
            mesh.nodes[i,:] = np.array(self.data.GetPoint(i),dtype=np.float64)
        # work on counting and assigning element types
        mesh.nel = self.data.GetNumberOfCells()
        num_hex8 = 0 # number of hex8
        num_tet4 = 0 # number of tet4
        # count total number of elements firstly
        for i in range(mesh.nel):
            if (self.data.GetCell(i).GetCellType()==12): # hexa8
                num_hex8 = num_hex8 + 1
            elif (self.data.GetCell(i).GetCellType()==10): #tet4
                num_tet4 = num_tet4 + 1
        # assign memory for mesh
        mesh.hex8 = np.zeros([num_hex8,8],dtype=np.uint64)
        mesh.tet4 = np.zeros([num_tet4,4],dtype=np.uint64)
        
        
        # read and assign element connectivities
        num_hex8 = 0 # number of hex8
        num_tet4 = 0 # number of tet4        
        for i in range(self.data.GetNumberOfCells()):
            if (self.data.GetCell(i).GetCellType()==12): # hexa8
                for j in range(self.data.GetCell(i).GetNumberOfPoints()):
                    mesh.hex8[num_hex8][j] = self.data.GetCell(num_hex8).GetPointId(j)
                num_hex8 = num_hex8 + 1
            elif (self.data.GetCell(i).GetCellType()==10): #tet4
                num_tet4 = num_tet4 + 1
                for j in range(self.data.GetCell(i).GetNumberOfPoints()):
                    mesh.tet4[num_tet4][j] = self.data.GetCell(num_tet4).GetPointId(j)
                num_tet4 = num_tet4 + 1
        # update nen       
        mesh.update_nel_nn()
        
        
        """
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
                hex = vtk.vtkHexahedron()
                for j in range(8):
                    hex.GetPointIds().SetId(j,mesh.hex8[i,j])
                self.uGrid.InsertNextCell(hex.GetCellType(),hex.GetPointIds())
    
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
    """    
    def __del__(self):
        class_name = self.__class__.__name__
        print (class_name, "destroyed")
        
