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
M_x = 100000.
M_y = 100000.
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

## Section 1
##for cut in range(np.shape(a)[0]-1):
##    x = a[cut]*np.cos(angle)
##    y = b[cut]*np.sin(angle)
##    I_xx = moi_hellipse(a[cut],b[cut])[0]
##    I_yy = moi_hellipse(a[cut],b[cut])[1]
##    sigma = M_x/I_xx*y+M_y/I_yy*x
##    for material in range(np.shape(mat.rho)[0]-1):
##    
##        t = sigma/mat.Fty[material]
##        t_maxt =  np.amax(t)
##        t_maxc =  -np.amin(t)
##        
##        if t_maxt > t_maxc:
##            t_material[cut][material] = t_maxt 
##        else:
##            t_material[cut][material] = t_maxc
##    if cut>1 and cut < np.shape(a)[0]-1:
##        t_const = np.maximum(t_material[cut],t_material[cut-1])
##
##for material in range(np.shape(mat.rho)[0]-1):
##    dm = mat.rho[material]*t_const[material]*dA*dl
##    mass[material] = sum(dm)
##    cost[material] = mass[material]*mat.Cost[material]

## To keep the cross-section symmetrical the number of booms is always kept at
## 4n. 
a = 0.15
b = 0.6
## Section 1 & 3
n_boom = 24
def moi_I(n_boom,a,b):
    angleboom = np.arange(0,360,360/n_boom)*np.pi/180
    x = a*np.cos(angleboom)
    y = b*np.sin(angleboom)
    I_xx = sum(y**2)
    I_yy = sum(x**2)
    return x,y,I_xx,I_yy

## Section 2
a_ellipse_top = 0.6
b_ellipse_top = 0.15
L_max         = 1.2
n_boom_top    = 9
n_boom_bottom = n_boom_top

angleboom_top = np.arange(0,180,180/n_boom_top)*np.pi/180
angleboom_top = np.append(angleboom_top,np.pi)

def f_x(x,a,b):
    f = (a**2*(np.sin(x)**2)+b**2*(np.cos(x)**2))**0.5
    return f
def arclength(angleboom_top,a,b):
    nodes = np.array([angleboom_top[0],angleboom_top[0]+ \
                 (angleboom_top[1]-angleboom_top[0])/3,\
                  angleboom_top[0]+2./3*(angleboom_top[1]-angleboom_top[0]),\
                  angleboom_top[1]])
    boom_spacing = np.around((angleboom_top[1]-angleboom_top[0])/8* \
                         (f_x(nodes[0],a,b)+ \
                        3*f_x(nodes[1],a,b)+ \
                        3*f_x(nodes[2],a,b)+ \
                        f_x(nodes[3],a,b)),1)
    return boom_spacing
dl = arclength(angleboom_top,a,b)
d = np.arange(0,L_max,dl)
d = L_max/2

x = a_ellipse_top*np.cos(angleboom_top)
print(x)
y = b_ellipse_top*np.sin(angleboom_top)+d
print(y)
while space < d/dl
    x  = x.append(x,x(np.shape[0]))
    print x
    y  = y.append(y,y([np.shape[0])-dl)
    
def boom_area(n_boom,a,b,M_x,M_y,x,y,I_xx,I_yy,materials):
    sigma_bending_B = M_x/I_xx*y+M_y/I_yy*x
    sigma_bending_B_max = np.amax(sigma_bending_B)
    sigma_bending_B_ind = np.argmax(sigma_bending_B)
    sigma_bending_B_min = -np.amas(sigma_bending_B)
    sigma_bending_B_ind = np.argamin(sigma_bending_B)
    
    if sigma_bending_B_max > sigma_bending_B_min:
        sigma_B = sigma_bending_B_max
    else:
        sigma_B = sigma_bending_B_min
    
    B = np.zeros(np.shape(materials)[0]-1)
    for material in range(np.shape(materials)[0]-1):
        B[material]= sigma_B/materials[material]*10**6
    return B

def skin_tickness(n_boom,a,b,S_x,S_y,T,x,y,I_xx,I_yy,materials):
    q_shear[0] = 0
    panel = 1
    while panel <= n_boom:
        q_shear[panel] = q_shear[panel-1]-S_x/I_yy*x-S_y/I_yy*y
        panel+=1
    A = np.pi*a*b
    q_torque = -T/(2*A)
    q = q_shear+q_torque
    q_max = np.amax(q)
    q_max_ind = np.argamax(q)
    q_min = -np.amin(q)
    q_min_ind = np.argmin(q)
    if q_max > q_min:
        q = q_max
    else:
        q = q_min
    t = np.zeros(np.shape(materials)[0]-1)
    for material in range(np.shape(materials)[0]-1):
        t[material]=q_max/materials[material]
    return t


   
    
    


