#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 22:38:55 2018

@author: jungjaeyong
"""


import os
import sys
file_dir = os.path.dirname(os.path.abspath("test_prepro.py"))
sys.path.append(file_dir+"/../src") # add relative path
#sys.path.append(os.getcwd()+"/../src") # add relative path
from mesh import *
from postp import *
from prepro import *

gmesh = mesh("global")
mypreprop = prepro()
mypreprop.read("test.vtk")
mypreprop.assign_geom(gmesh)
print("preprocessing done")
mypos = postp(gmesh)
print("postprocessing done")
mypos.write("test2.vtk")
print("writing done")

