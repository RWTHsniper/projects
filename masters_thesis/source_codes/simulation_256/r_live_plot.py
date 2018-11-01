import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import operator
import numpy as np


fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

def animate(i):
    data = [[],[]]
    with open(filename) as f:
	for line in f:
	    dum_data = line.split()
	    for i in range(2):
		data[i].append(float(dum_data[i]))

    maxindex, maxim = max(enumerate(data[1]), key=operator.itemgetter(1))
    minim = min(data[1])

    ax1.clear()
    ax1.plot(data[0],data[1],color="pink",label="Stress-Strain curve",linewidth=5)
    plt.legend(loc='upper right')
    plt.title(title)
    plt.xlabel("Strain",fontsize=15)
    plt.ylabel("Stress",fontsize=15)
    plt.text(data[0][maxindex], data[1][maxindex], r'Peak ('+str(data[0][maxindex])+", "+str(data[1][maxindex])+")")
    len_data = [(data[0][-1]-data[0][0]),(maxim-data[1][0])]
    my_axis = [(data[0][0]-fct*len_data[0]), (data[0][-1]+fct*len_data[0]), (data[1][0]-fct*len_data[1]), (maxim+fct*len_data[1])]	# xmin, xmax, ymin, ymax
    ax = fig.gca()
    d_grid = [(my_axis[1]-my_axis[0])/(num_grids[0]+1),(my_axis[3]-my_axis[2])/(num_grids[1]+1)]	# grid space
    ax.set_xticks(np.arange(my_axis[0], my_axis[1],d_grid[0]))
    ax.set_yticks(np.arange(my_axis[2], my_axis[3],d_grid[1]))
    plt.grid()
    plt.axis(my_axis)
    
def init(filename,title,num_grids,fct):
    0  
# Set file name to read
#    filename = filename  
#    title = title

    
filename = "r_z_stress_strain.txt"
title = "restart_Composite_NH_Stress-Strain(exponential damage function)"
num_grids=[6,5]	# number of grids in each dimension
fct=0.1	# Factor of offset for outerskirt
    
#ani = animation.FuncAnimation(fig, animate, init_func=init(filename), interval=1000)
#init(filename,title)
ani = animation.FuncAnimation(fig, animate, init_func=init(filename,title,num_grids,fct), interval=1000)
plt.show() 
fig.savefig(title+'.jpg')


