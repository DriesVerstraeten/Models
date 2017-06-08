# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:54:31 2017

@author: 123
"""

#Calculations for the FUSELAGE MODEL
# Input: Semi-major axis of the thin walled ellipse (a)
#        Semi-minor axis of the thin walled ellipse (b)
#        Thickness of the ellipse (t)
# Output: MOI x-axis
#         MOI y-axis
# The coordinate system has its orgin always at the end of a section, with the
# x-axis pointing in the positive x direction. 
###############################################################################
import numpy as np 
import Material_properties as mat

a = 2
b = 1 
x = 1
y = 2
M_x = 10000
M_y = 10000
    
I_xx = []
I_yy = []
    
I_xx.append(np.pi / 4 * b**2 *(3*a+b)) # Without thickness
I_yy.append(np.pi / 4 * a**2 *(3*b+a)) # Without thickness
p = np.shape(rho) 
t = [];
#for i in range(np.shape(rho)[0]-1) 
#    t.append((M_x/(I_xx)*y + M_y/I_yy*x)/Ftu[i])
    