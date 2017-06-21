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
        I_xx = moi_circle(r[i])
        I_yy = I_xx
        sigma_t = M_x[i]/I_xx*y+M_y[i]/I_yy*x 
        sigma_t_max = np.amax(sigma_t)
        sigma_t_min = -np.amin(sigma_t)
        if sigma_t_max > max_sigma_t:
            max_sigma_t = sigma_t_max
        if sigma_t_min > min_sigma_t:
            min_sigma_t = sigma_t_min
    return max_sigma_t, min_sigma_t

def shear_i(r,S_x,S_y,T):
    max_tau_t = 0
    min_tau_t = 0
    for i in range(np.shape(r)[0]):
        x = section_i(r[i])[0]
        y = section_i(r[i])[1]
        I_xx = moi_circle(r[i])
        I_yy = I_xx
        q_b  = r[i]*(-S_x[i]/I_yy*y+S_y[i]/I_xx*(x-r[i]))
        q_s0 = T[i]/(2*np.pi*r[i])
        q    = q_b+q_s0
        q_max = np.amax(q)
        q_min = -np.amin(q)
        if q_max > max_tau_t:
            max_tau_t = q_max
        if q_min > min_tau_t:
            min_tau_t = q_min      
    if max_tau_t > min_tau_t:
        t = max_tau_t
    else:
        t = min_tau_t
    return t

def ip(r,dp,sigma_hoop):
    t = dp*r/sigma_hoop
    t = np.amax(t)
    return t

def thickness(r,M_x,M_y,S_x,S_y,T,sigma_t,sigma_c,tau,t_c_min,dp,sigma_hoop,rho_f,rho_c,E_x,E_c):
    t_bt = bending_i(M_x,M_y,r)[0]/sigma_t
    print(t_bt)
    t_bc = bending_i(M_x,M_y,r)[1]/sigma_c
    t_s  = shear_i(r,S_x,S_y,T)/tau
    print(t_s)
    t_h  = ip(r,dp,sigma_hoop)
    t    = np.array([t_bt, t_bc, t_s, t_h])
    t_eq = np.round(np.amax(t),3)
    print(t_eq)
    #t_f  = t_eq/2*E_x/(71.7*10**9)
    #t_f   = np.round(t_eq/2*23.1/(26),4)
    #delta = 36*t_f**2-4*(4*t_f**2-t_eq**2)*3
    #delta = 36*t_f**2-4*(4*t_f**2-t_eq**3*71.7*10**9/(2*t_f*E_x))*3
    #h_c = (-6*t_f + delta**(0.5))/(2*3)
    #h_c = t_eq-2*t_f
    h_c = t_c_min
    t_f = np.round((t_eq*rho_f-h_c*rho_c)/(2*rho_f),3)
    #t_f = (t_eq-h_c)/2
    return t_f, h_c

def mass(r,M_x,M_y,S_x,S_y,T,sigma_t,sigma_c,tau,t_c_min,dp,sigma_hoop,rho_f,rho_c,E_x,E_c):
        dm = (2*thickness(r,M_x,M_y,S_x,S_y,T,sigma_t,sigma_c,tau, \
                          t_c_min,dp,sigma_hoop,rho_f,rho_c,E_x,E_c)[0] \
         *rho_f+thickness(r,M_x,M_y,S_x,S_y,T,sigma_t,sigma_c, \
                          tau,t_c_min,dp,sigma_hoop,rho_f,rho_c,E_x,E_c)[1]*\
         rho_c)*2*np.pi*r*0.01
        return dm




