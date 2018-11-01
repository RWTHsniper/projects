
# I am making live graph
import numpy as np
import pylab as plt
import operator
import os

#import matplotlib.animation as animation
#from matplotlib import style

"""
Here, I will write down textfile reader for real number data.
"""

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

# Parameters for graph
num_grids=[6,5]	# number of grids in each dimension
fct=0.1	# Factor of offset for outerskirt

# Title should be defined
title = "restart_Composite_NH_Stress-Strain(exponential damage function)"
# Import data
a = "r_z_stress_strain.txt"
reader = file_reader(a,2)
data = reader.return_data()
#maxim, maxindex = find_max(sum_data[1])
#index, value = max(enumerate(my_list), key=operator.itemgetter(1))
# Find maximum
maxindex, maxim = max(enumerate(data[:,1]), key=operator.itemgetter(1))
minim = min(data[:,1])
fig = plt.figure(1)
plt.plot(data[:,0],data[:,1],color="pink",label="Stress-Strain curve",linewidth=5)
plt.legend(loc='upper right')
#fig.suptitle("Stress-Strain")
plt.title(title)
plt.xlabel("Strain",fontsize=15)
plt.ylabel("Stress",fontsize=15)
plt.text(data[maxindex][0], data[maxindex][1], r'Peak ('+str(data[maxindex][0])+", "+str(data[maxindex][1])+")")

len_data = [(data[-1][0]-data[0][0]),(maxim-data[0][1])]
my_axis = [(data[0][0]-fct*len_data[0]), (data[-1][0]+fct*len_data[0]), (data[0][1]-fct*len_data[1]), (maxim+fct*len_data[1])]	# xmin, xmax, ymin, ymax

## set grid size
ax = fig.gca()
#ax.set_xticks(np.arange(data[0][0], data[0][-1], (data[0][-1]-data[0][0])/num_grids))
#ax.set_yticks(np.arange(minim, maxim, (maxim-minim)/num_grids))
d_grid = [(my_axis[1]-my_axis[0])/(num_grids[0]+1),(my_axis[3]-my_axis[2])/(num_grids[1]+1)]
# grid space
ax.set_xticks(np.arange(my_axis[0], my_axis[1],d_grid[0]))
ax.set_yticks(np.arange(my_axis[2], my_axis[3],d_grid[1]))
plt.grid()
plt.axis(my_axis)

#plt.axis([0, 1, 0, 1])
plt.show()
fig.savefig(title+'.jpg')






