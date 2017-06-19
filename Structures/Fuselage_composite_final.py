# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:18:06 2017

@author: Claudia
"""
## Composite Analysis:
import numpy as np
#import Fuselage_parameters as fp

def moi_circle(r):
    A = 2*np.pi*r
    I = np.pi*r**3
    return I

def arclength(angleboom,R):
    nodes = angleboom[1]
    spacing = 2*np.pi*R*(nodes)/(2*np.pi)
    return spacing

def section_i(r):
    angle = np.arange(360)*np.pi/180
    x = r*np.cos(angle)
    y = r*np.sin(angle)
    return x,y
    
def bending_i(M_x,M_y,r):
    max_sigma_t = 0
    min_sigma_t = 0 
    for i in range(np.shape(r)[0]):
        x = section_i(r[i])[0]
        y = section_i(r[i])[1]
        I_xx = moi_circle(r[i])[0]
        I_yy = I_xx
        sigma_t = M_x/I_xx*y+M_y/I_yy*x 
        sigma_t_max = np.amax(sigma_t)
        sigma_t_min = -np.amax(sigma_t)
        if sigma_t_max > max_sigma_t:
            max_sigma_t = sigma_t_max
        if sigma_t_min > min_sigma_t:
            min_sigma_t = sigma_t_min
    return max_sigma_t,min_sigma_t

def shear_i(r,S_x,S_y,T):
    max_tau_t = 0
    for i in range(np.shape(r)[0]):
        x = section_i(r[i])[0]
        y = section_i(r[i])[1]
        I_xx = moi_circle(r[i])[0]
        I_yy = I_xx
        q_b  = r[i]*(-S_x/I_yy*y+S_y/I_xx*(x-r[i]))
        q_s0 = T/(2*np.pi*r[i])
        q    = q_b+q_s0
        q_max = np.amax(q)
        if q_max > max_tau_t:
            max_tau_t = q_max
    return max_tau_t 

def ip(r,dp,sigma_hoop):
    t = dp*r/sigma_hoop
    return t

def thickness(r,M_x,M_y,S_x,S_y,T,sigma_t,tau,t_c_min,dp,sigma_hoop):
    t_bt = bending_i(M_x,M_y,r)[0]/sigma_t
    t_bc = bending_i(M_x,M_y,r)[1]/sigma_c
    t_s = shear_i(r,S_x,S_y,T)/tau
    t_h = ip(r,dp,sigma_hoop)
    t   = np.array([t_bt,t_bc,t_s,t_h])
    t_eq = np.amax(t)
    t_f  = t_eq/2
    delta = 36*t_f**2-4*(4*t_f**2-t_eq)*3
    h_c = (-6*t_f + delta**(0.5))/(2*3)
    if h_c < 0:
        h_c = t_c_min
    return t_f,h_c
    
    
##for material in range(np.shape(mat.rho)[0]-1):
##    dm = mat.rho[material]*t_const[material]*dA*dl
##    mass[material] = sum(dm)
##    cost[material] = mass[material]*mat.Cost[material]
##
    



