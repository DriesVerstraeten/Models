# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:57:36 2017

@author: 123
"""

#Calculations for the HORIZ VERT TAILS MODEL

import Parameters as p
import Material_properties as mat
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import integrate

#WINGBOX DIMENSIONS
box_t = 0.0005 #m, thickness. same everywhere for now
box_b = 0.25*p.cr_ht #m, width at rootapply scaling for other points. scales linearly with taper
box_h = 0.1*p.cr_ht #m, height at root. apply scaling for other points. scales linearly with taper

def y_tail(dy,i):
    y_tail = np.linspace(0,p.b_ht/2,dy)
    return y_tail[i]

def y_tail_inv(dy,i):
    y_tail_inv = np.linspace(p.b_ht/2,0,dy)
    return y_tail_inv[i]

def c_tail(dy,i):
    c_tail = np.linspace(p.ct_ht,p.cr_ht,dy)
    return c_tail[i]
    
def rho_tail(dy,i):
    mid_rho = p.W_ht*2/((p.cr_ht+p.ct_ht)*p.b_ht)
    root_rho = mid_rho / 0.9
    tip_rho = root_rho * 0.8
    rho = np.linspace(tip_rho, root_rho, dy)
    return rho[i]

def area(dy,i):
    areas = p.b_ht/(2*dy) * c_tail(dy,i)
    return areas

def slope_c(dy):
    slope = (p.ct_ht+p.cr_ht)/dy
    return slope

def tail_I_xx(dy,i):
    I_xx = np.zeros(dy)
    I_xx = (p.cr_ht - y_tail_inv(dy,i) * slope_c(dy))**4 * (box_h)**3 * box_b/12 - (p.cr_ht - y_tail_inv(dy,i) * slope_c(dy))**4 * (box_h - box_t)**3 * (box_b - 0.0005)/12
    return I_xx
    #9G ULTIMATE STATIC LOAD HORIZONTAL TAIL

#Distribution at every point
def tail_force(dy,g,i):
    force = g*p.g*rho_tail(dy,i)*area(dy,i)
    return force

def tail_distr(dy,g,i):
    distr = g*p.g*rho_tail(dy,i)*area(dy,i)/(p.b_ht/2000)
    return distr 

def tail_shear(dy,g,j):
    shear = np.zeros(j+1)
    shear[0] = tail_force(dy,g,0)
    for i in range (1,j+1):
        shear[i] = shear[i-1] + tail_force(dy,g,i)
    return shear[j]

def tail_moment(dy,g,j):
    moment = np.zeros(j+1)
    moment[0] = tail_force(dy,g,0) * y_tail(dy,0)
    for i in range (1,j+1):
        moment[i] = moment[i-1] + tail_force(dy,g,i) * y_tail(dy,i)
    return moment[j]
        
def tail_plots(dy,g):
    distr = np.zeros(dy)
    y = np.zeros(dy)
    shear = np.zeros(dy)
    moment = np.zeros(dy)
    for i in range (0,dy):
        distr[i] = tail_distr(dy,g,i)
    for i in range (0,dy):
        y[i] = y_tail(dy,i)
    for i in range (0,dy):
        shear[i] = tail_shear(dy,g,i)
    for i in range (0,dy):
        moment[i] = tail_moment(dy,g,i)
    plt.figure(figsize=(19,5))
    plt.subplot(131)
    plt.plot(y, distr)
    plt.ylabel('Load, N')
    plt.xlabel('Location, m')
    plt.subplot(132)
    plt.plot(y, shear)
    plt.ylabel('Shear, N')
    plt.xlabel('Location, m')
    plt.subplot(133)
    plt.plot(y, moment)
    plt.ylabel('Moment, Nm')
    plt.xlabel('Location, m')
    plt.show()
    return

def tail_shear_stress(dy,dx,dz,i,g):
    x_12 = np.linspace(0,((p.cr_ht - y_tail_inv(dy,i) * slope_c(dy)) * box_b/2),dx)
    z_23 = np.linspace(((p.cr_ht - y_tail_inv(dy,i) * slope_c(dy)) * box_h)/2, -1*(p.cr_ht-y_tail_inv(dy,i)*slope_c(dy))*0.1*p.cr_ht/2, dz)
    x_34 = np.linspace((p.cr_ht - y_tail_inv(dy,i) * slope_c(dy)) * box_b/2,0,dx)
    q_12 = np.zeros(dx)
    q_23 = np.zeros(dy)
    q_34 = np.zeros(dx)
    for j in range (0,dx):
        q_12[j] = -tail_shear(dy,g,i) * box_t * z_23[0] * x_12[j] / (2*tail_I_xx(dy,i))
    for j in range (0,dz):
        q_23[j] = q_12[dx-1] + box_t * tail_shear(dy,g,i) * x_12[dx-1] * ((z_23[0])**2 - z_23[j]**2)/(4*tail_I_xx(dy,i))
    for j in range (0,dx):
        q_34[j] = q_23[dz-1] - tail_shear(dy,g,i) * box_t * z_23[dz-1] * x_12[j] / (2*tail_I_xx(dy,i))
    tau_12 = np.zeros(dx)
    tau_23 = np.zeros(dx)
    tau_34 = np.zeros(dx)
    for j in range (0,dx):
        tau_12[j] = q_12[j] / box_t
    for j in range (0,dx):
        tau_23[j] = q_23[j] / box_t
    for j in range (0,dx):
        tau_34[j] = q_34[j] / box_t
    plt.figure(figsize=(19,5))
    plt.subplot(231)
    plt.plot(x_12, q_12)
    plt.ylabel('Shear, N/m')
    plt.xlabel('Location, m')
    plt.subplot(232)
    plt.plot(q_23, z_23)
    plt.ylabel('Location, m')
    plt.xlabel('Shear, N/m')
    plt.subplot(233)
    plt.plot(x_34, q_34)
    plt.ylabel('Shear, N/m')
    plt.xlabel('Location, m')
    plt.subplot(234)
    plt.plot(x_12, tau_12)
    plt.ylabel('Shear stress, N/m^2')
    plt.xlabel('Location, m')
    plt.subplot(235)
    plt.plot(z_23, tau_23)
    plt.ylabel('Location, m')
    plt.xlabel('Shear, N/m')
    plt.subplot(236)
    plt.plot(x_34, tau_34)
    plt.ylabel('Shear, N/m')
    plt.xlabel('Location, m')
    plt.show()   
    return













#Deflection
#a = 0.5 * forces_9g[999]
#b = 0.25 * (distr_9g[999] - distr_9g[998])/p.b_ht/2000
#c = p.cr_ht
#d = slope_c
#e = mat.E[0] * (0.1*p.cr_ht)**3 * 0.25*p.cr_ht/12 - (0.1*p.cr_ht - 0.0005)**3 * (0.25*p.cr_ht/12 - 0.0005)
#v_9g = (-(b*d*y_tail+a*d-4*b*c)*np.log(abs(d*y_tail-c)))/(e*d*d*d*d*d) + (c*(3*(2*a*d-3*b*c)*d*y_tail-(5*a*d-8*b*c)*c))/(6*d*d*d*d*d*e*(d*y_tail-c)*(d*y_tail-c)) + (b*y_tail)/(d*d*d*d*e) + y_tail*(b*np.log(abs(-c)))/(e*d*d*d*d) - y_tail*(2*a*d-11*b*c)*c*c/(6*d*d*d*d*c*c*c*e) + (a*d-4*b*c)*np.log(abs(-c))/(d*d*d*d*d*e) + ((5*a*d-8*b*c)*c)/(6*d*d*d*d*d*c*c*e)       
#plt.subplot(144)
#plt.plot(y_tail, v_9g)
#plt.show()       
#Deflection
#deflec_45g = forces_45g * (p.b_ht/4000)**2 * (3 * p.b_ht/2000-p.b_ht/4000) /(6 * mat.E[0] * I_xx)
#for i in range (2,dy+1):
#    deflec_45g[dy-i] = deflec_45g[dy-i] + deflec_45g[dy-i+1]
#plt.subplot(144)
#plt.plot(y_tail, deflec_45g)
#plt.show()
