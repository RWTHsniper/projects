# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 14:55:53 2017

@author: jjung

Modified on Fri Nov 4
Ability to write Eta1 for C/C++ is added.

"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import random as rand
#from pyevtk.hl import gridToVTK

class grid:
   'Indexing of grid is 0-base'

   def __init__(self, L,N):
       self.L = L
       self.N = N
       self.Eta1 = np.zeros((N[0],N[1],N[2]),dtype=np.int,order='C')

       dx = [0,0,0]
       for i in range(len(L)):
           dx[i] = L[i]/(float(N[i]))
       self.dx = dx
       self.midpoint = 0.5*np.array([dx[0]*(N[0]-1.0),dx[1]*(N[1]-1.0),dx[2]*(N[2]-1.0)])
       points = np.zeros((3,N[0],N[1],N[2]),dtype=np.float64,order='C')
       for k in range(N[2]):
           for j in range(N[1]):
               for i in range(N[0]):
                   points[0,i,j,k] = i*dx[0]
                   points[1,i,j,k] = j*dx[1]
                   points[2,i,j,k] = k*dx[2]
       self.x = points
       self.npoints = N[0]*N[1]*N[2]
       self.cellpoints = (N[0]-1)*(N[1]-1)*(N[2]-1)

      # calculate grid point coordinate according to grid index
   def get_position(self,indx):
       position = [self.x[0,indx[0],indx[1],indx[2]],self.x[1,indx[0],indx[1],indx[2]],self.x[2,indx[0],indx[1],indx[2]]]
       return position
      # Assign material property on each grid point
   def assign_prop(self,geom_obj):
       for k in range(self.N[2]):
           for j in range(self.N[1]):
               for i in range(self.N[0]):
                   pnt = self.get_position([i,j,k])
                   chk = geom_obj.check_in_out(pnt)
                   if ((self.Eta1[i,j,k] != 0)and(chk == 0)):
                       # Eta1 is assigned and chk is just zero
                       continue
                   if (( self.Eta1[i,j,k] != chk )and(self.Eta1[i,j,k] != 0)and(chk != 0)):
                       # Material assigned already but chk is nonzero
                       print "Eta1 is overlapping at ", i, " ", j, " ", k
                   self.Eta1[i,j,k] = chk
      # Assign material property on each grid point according to periodic B.C.
   def periodic_assign_prop(self,geom_obj):
       flag_good_2_assign = 1
       for k in range(self.N[2]):
           if (flag_good_2_assign == 0):
               break;
           for j in range(self.N[1]):
               if (flag_good_2_assign == 0):
                  break;
               for i in range(self.N[0]):
                   # if material is already assigned, do not go through assignment
                   if (flag_good_2_assign == 0):
                       break;
		   if ((self.Eta1[i,j,k] != 0)):
                       chk = geom_obj.periodic_check_in_out(self.get_position([i,j,k]),self.midpoint)
#                   if (( self.Eta1[i,j,k] != chk )and(self.Eta1[i,j,k] != 0)and(chk != 0)):
                       if ((chk != 0)):
                       # Material assigned already but chk is nonzero
                           print "Eta1 is overlapping at ", i, " ", j, " ", k
                           flag_good_2_assign = 0
		           return
       if (flag_good_2_assign == 1):
           for k in range(self.N[2]):
               for j in range(self.N[1]):
                   for i in range(self.N[0]):
		       if (self.Eta1[i,j,k] == 0):
                           chk = geom_obj.periodic_check_in_out(self.get_position([i,j,k]),self.midpoint)
                           if ((self.Eta1[i,j,k] == 0)and(chk != 0)):
                               self.Eta1[i,j,k] = chk




   def plot_3dgrid(self,title,save_fig,show_fig):
       colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
       markers = ['.',',','o','^']
       markers_size = [70,50,50,50,50]
       fig = plt.figure()
       ax = Axes3D(fig)
       if ((save_fig + show_fig) > 0):
           for k in range(self.N[2]):
               for j in range(self.N[1]):
                   for i in range(self.N[0]):
                       pnt = self.get_position([i,j,k])
                       if ((self.Eta1[i,j,k] != 0)):
		           ax.scatter(pnt[0], pnt[1], pnt[2], s=markers_size[self.Eta1[i,j,k]], c=colors[self.Eta1[i][j][k]], marker=markers[self.Eta1[i][j][k]])
                       elif (((i == 0)or(i == self.N[0]-1))and((j == 0)or(j == self.N[1]-1))and((k == 0)or(k == self.N[2]-1))):
                           ax.scatter(pnt[0], pnt[1], pnt[2], s=markers_size[self.Eta1[i,j,k]], c=colors[self.Eta1[i][j][k]], marker=markers[self.Eta1[i][j][k]])
           ax.set_xlabel('X Label')
           ax.set_ylabel('Y Label')
           ax.set_zlabel('Z Label')
       if (save_fig == 1):
           fig.savefig(title+".png", bbox_inches='tight')
       if (show_fig == 1):
           plt.show()

      # Destroyer
   def __del__(self):
       class_name = self.__class__.__name__
       print class_name, "destroyed"


class cylinder:
   'Common base class for all employees'

   def __init__(self,L, center,direction, radius, height, prop):
      self.grid_len = L
      self.center = center
      direction = direction / np.linalg.norm(direction)
      self.direction = np.array(direction);
      self.direction /= np.linalg.norm(direction) # normalized
      print "direction",self.direction
      self.radius = radius
      self.height = height
      self.mat_prop = prop
      # Check if a grid point is inside or outside
      # return 1: It is inside cylinder
   def check_in_out(self,point):
       tol = 1.0e-3
       flag = 0
       r = [a - b for a, b in zip(point, self.center)]
       d1 = np.dot(self.direction,r)
       if (d1 < 0.0):
           direct = -self.direction
           d1 = -d1
       else:
           direct = self.direction
       if (d1 <= (self.height/2.0)*(1.0+tol)):
           r2 = [a - d1*b for a, b in zip(r, direct)]
           d2 = np.linalg.norm(r2)
           if (d2 < self.radius):
               flag = self.mat_prop
       return flag
   # check in-out based on periodic B.C.
   def periodic_check_in_out(self,point,midpoint):
       tol = 1.0e-7
       shifter = np.zeros(3,dtype=np.float64)
       shift_cent = np.zeros(3,dtype=np.float64)
# check in-out by shifting center of the object
# x-direction
       if (self.center[0] < midpoint[0]):
           shifter[0] = self.grid_len[0]
       else:
           shifter[0] = -self.grid_len[0]
# y-direction
       if (self.center[1] < midpoint[1]):
           shifter[1] = self.grid_len[1]
       else:
           shifter[1] = -self.grid_len[1]
# z-direction
#       if (self.center[2] < midpoint[2]):
#           shifter[2] = self.grid_len[2]
#       else:
#           shifter[2] = -self.grid_len[2]
       for i in range(2):
           for j in range(2):
#               for k in range(2):
#		   print "i,j,k",i,j,k
               shift_cent = self.center + np.array([shifter[0]*i,shifter[1]*j,0.0])
#As infinite length of fiber is applied. There should be no shift consideration in z-direction"
#		   print"direction",self.direction
#                   print "shift_center",shift_cent
               r = [a - b for a, b in zip(point, shift_cent)]
               d1 = np.dot(self.direction,r)
               if (d1 < 0.0):
                   direct = -self.direction
                   d1 = -d1
               else:
                   direct = self.direction
               if (d1 <= (self.height/2.0)):
                   r2 = [a - d1*b for a, b in zip(r, direct)]
                   d2 = np.linalg.norm(r2)
                   if (d2 <= self.radius):
		       return self.mat_prop

       return 0



      # Destroyer
   def __del__(self):
       class_name = self.__class__.__name__
       print class_name, "destroyed"

class sphere:
   'Common base class for all employees'
   def __init__(self, center, radius, prop):
      self.center = center
      self.radius = radius
      self.mat_prop = prop
      # Check if a grid point is inside or outside
      # return 1: It is inside cylinder
   def check_in_out(self,point):
       flag = 0
       r = [a - b for a, b in zip(point, self.center)]
       d1 = np.linalg.norm(r)
       if (d1 <= self.radius):
           flag = self.mat_prop
       return flag
      # Destroyer
   def __del__(self):
       class_name = self.__class__.__name__
       print class_name, "destroyed"


def write_Eta1(outpfile, Eta1, N,wtype):
    file = open(outpfile,"w")
    if (wtype =="F"):
        for k in range(N[2]):
            for j in range(N[1]):
                for i in range(N[0]):
                    file.write(str(Eta1[i][j][k])+"\n")
    elif(wtype == "C"):
        for i in range(N[0]):
            for j in range(N[1]):
                for k in range(N[2]):
                    file.write(str(Eta1[i][j][k])+"\n")
    else:
        print "Type of Eta1 file should be defined."
    file.close()
"""
How to use this script.
1. define string of title
>>  title1 = "Eta1_123"
2. define grid object using length of grid and number of grid points.
>> L = [Lx,Ly,Lz]; N = [Nx,Ny,Nz]
>> grid1 = grid(L,N)
3. define object of geometry. For current version, only cyliner object is done and good to go.
  define cylinder object using (Length array of grid, direction vector, radius, length of cylinder, assigned Eta1 value)
>>   center = [cx,cy,cz]
>>   direction_vector = [0,0,1.0]
>>   cylin1 = cylinder(L,center,direction_vector,radius,Lz,(1))
4. Assign cylinder object on grid object. A cylinder object on a boundary of RVE will have periodic continuity over the boundary.
>>   grid1.periodic_assign_prop(cylin1)
>>   del cylin1
5. Write_Eta1 file.
>>  write_Eta1(title1, grid1.Eta1, grid1.N)
"""




# Please check this input section before running
title1 = "Eta1_128"
title2 = "Eta1_256"
title3 = "Eta1_64"
save_fig = 0; show_fig = 0
Lx = 100.0e-3;Ly = 100.0e-3;Lz = 100.0e-3
Nx = 128; Ny = 128; Nz = 1
num_inclusions = 60;
L = [Lx,Ly,Lz]; N = [Nx,Ny,Nz]
############################################
grid1 = grid(L,N)
N=[256,256,1]
grid2 = grid(L,N)
N=[64,64,1]
grid3 = grid(L,N)
for i in range(num_inclusions):
  rx=rand.randrange(0,Nx)
  ry=rand.randrange(0,Ny)
  rand_radi = rand.random()/3
  rx = Lx*rx/Nx;   ry = Ly*ry/Ny
  print "Process", str(i+1)+"/"+str(num_inclusions)
  print "rx,ry,rand_radi", rx, ry,rand_radi
  cylin1 = cylinder(L,[rx,ry,L[2]/2.0],[0,0,1.0],(1+rand_radi)*Lx/16.0,3.0*Lz,(1))
  grid1.periodic_assign_prop(cylin1)
  grid2.periodic_assign_prop(cylin1)
  grid3.periodic_assign_prop(cylin1)
  del cylin1

"""
sph1 = sphere([Lx/2.0,Ly/2.0,Lz/2.0],Lx/4.0,2)
cylin1 = cylinder(L,[0,0,Lz/2.0],[0.0,0.0,1.0],Lx/12.0,Ly,1)
grid1.periodic_assign_prop(cylin1)
del cylin1
 """

#grid1.plot_3dgrid(title,1,1)
write_Eta1(title1, grid1.Eta1, grid1.N,"C")
write_Eta1(title2, grid2.Eta1, grid2.N,"C")
write_Eta1(title3, grid3.Eta1, grid3.N,"C")
#grid1.plot_3dgrid(title,save_fig,show_fig)
#gridToVTK("./"+title1, grid1.x[0], grid1.x[1], grid1.x[2], pointData = {"Eta1" : grid1.Eta1})
#gridToVTK("./"+title2, grid2.x[0], grid2.x[1], grid2.x[2], pointData = {"Eta1" : grid2.Eta1})
#rint "saving VTK file is done"
