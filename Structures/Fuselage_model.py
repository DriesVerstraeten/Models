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

def moi_hellipse(a,b):
    I_xx = np.pi / 4 * b**2 *(3*a+b)
    I_yy = np.pi / 4 * a**2 *(3*b+a)
    return I_xx,I_yy

a1 = 0.4
b1 = 0.3
M_x = 10000000.
M_y = 10000000.
L_1 = 1.2+0.15
a2 = 3.
b2 = 2.

L= 1.20
dl = 0.01

a = np.linspace(a1,a2,L/dl)
b = np.linspace(b1,b2,L/dl)
dA = a+b
angle = np.arange(361)*np.pi/180
t_material = np.zeros((np.shape(a)[0]-1,np.shape(mat.rho)[0]-1))
mass = np.zeros(np.shape(mat.rho)[0]-1)
cost = np.zeros(np.shape(mat.rho)[0]-1)
for cut in range(np.shape(a)[0]-1):
    x = a[cut]*np.cos(angle)
    y = b[cut]*np.sin(angle)
    I_xx = moi_hellipse(a[cut],b[cut])[0]
    I_yy = moi_hellipse(a[cut],b[cut])[1]
    sigma = M_x/I_xx*y+M_y/I_yy*x
    for material in range(np.shape(mat.rho)[0]-1):
    
        t = sigma/mat.Fty[material]
        t_maxt =  np.amax(t)
        t_maxc =  -np.amin(t)
        
        if t_maxt > t_maxc:
            t_material[cut][material] = t_maxt 
        else:
            t_material[cut][material] = t_maxc
    if cut>1 and cut < np.shape(a)[0]-1:
        t_const = np.maximum(t_material[cut],t_material[cut-1])

for material in range(np.shape(mat.rho)[0]-1):
    dm = mat.rho[material]*t_const[material]*dA*dl
    mass[material] = sum(dm)
    cost[material] = mass[material]*mat.Cost[material]

