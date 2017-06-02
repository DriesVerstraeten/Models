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

dy = 1000 #number of "pieces"
y_tail = np.linspace(0,p.b_ht/2,dy)
y_tail_inv = np.linspace(p.b_ht/2,0,dy)
c_tail = np.linspace(p.ct_ht,p.cr_ht,dy)
mid_rho = p.W_ht*2/((p.cr_ht+p.ct_ht)*p.b_ht)
root_rho = mid_rho / 0.9
tip_rho = root_rho * 0.8
rho_tail = np.linspace(tip_rho, root_rho, dy)
areas = p.b_ht/(2*dy) * c_tail
slope_c = (p.ct_ht+p.cr_ht)/dy
I_xx = np.zeros(dy)
I_xx = (p.cr_ht - y_tail_inv * slope_c)**4 * (0.1*p.cr_ht)**3 * 0.25*p.cr_ht/12 - (p.cr_ht - y_tail_inv * slope_c)**4 * (0.1*p.cr_ht - 0.0005)**3 * (0.25*p.cr_ht/12 - 0.0005)
#9G ULTIMATE STATIC LOAD HORIZONTAL TAIL

#Distribution at every point
forces_9g = np.zeros(1000)
distr_9g = np.zeros(1000)
for i in range (0,len(y_tail)):
   forces_9g[i] = 9*p.g*rho_tail[i]*areas[i]
for i in range (0,len(y_tail)):
   distr_9g[i] = 9*p.g*rho_tail[i]*areas[i]/(p.b_ht/2000)
plt.figure(figsize=(19,5))
plt.subplot(141)
plt.plot(y_tail, distr_9g)
plt.ylabel('Load, N')
plt.xlabel('Location, m')

#Shear distribution
shear_9g = list(forces_9g)
for i in range (1,len(y_tail)):
   shear_9g[i] = shear_9g[i-1] + forces_9g[i]
plt.subplot(142)
plt.plot(y_tail, shear_9g)
plt.ylabel('Shear, N')
plt.xlabel('Location, m')

#Moment at every point
moment_9g = forces_9g * y_tail
for i in range (1,len(y_tail)):
   moment_9g[i] = moment_9g[i-1]+moment_9g[i]
plt.subplot(143)
plt.plot(y_tail, moment_9g)
plt.ylabel('Moment, Nm')
plt.xlabel('Location, m')
plt.show()

#Deflection
a = 0.5 * forces_9g[999]
b = 0.25 * (distr_9g[999] - distr_9g[998])/p.b_ht/2000
c = p.cr_ht
d = slope_c
e = mat.E[0] * (0.1*p.cr_ht)**3 * 0.25*p.cr_ht/12 - (0.1*p.cr_ht - 0.0005)**3 * (0.25*p.cr_ht/12 - 0.0005)
v_9g = (-(b*d*y_tail+a*d-4*b*c)*np.log(abs(d*y_tail-c)))/(e*d**5) + (c*(3*(2*a*d-3*b*c)*d*y_tail-(5*a*d-8*b*c)*c))/(6*d**5*e*(d*y_tail-c)**2) + (b*y_tail)/(d**4*e) + y_tail*(b*np.log(-c))/d**4 - y_tail*(2*a*d-8*b*c)*c/(6*d**4*c**3) + (a*d-4*b*c)*/d**5 + ((5*a*d-8*b*c)*c)/(6*d**5*c**2)       
plt.subplot(144)
plt.plot(y_tail, v_9g)
plt.show()       
#-4.5G ULTIMATE STATIC LOAD HORIZONTAL TAIL

#Distribution at every point
forces_45g = np.zeros(1000)
distr_45g = np.zeros(1000)
for i in range (0,len(y_tail)):
   forces_45g[i] = -4.5*p.g*rho_tail[i]*areas[i]
for i in range (0,len(y_tail)):
   distr_45g[i] = -4.5*p.g*rho_tail[i]*areas[i]/(p.b_ht/2000)
plt.figure(figsize=(19,5))
plt.subplot(141)
plt.plot(y_tail, distr_45g)
plt.ylabel('Load, N/m')
plt.xlabel('Location, m')

#Shear distribution
shear_45g = list(forces_45g)
for i in range (1,len(y_tail)):
   shear_45g[i] = shear_45g[i-1] + forces_45g[i]
plt.subplot(142)
plt.plot(y_tail, shear_45g)
plt.ylabel('Shear, N')
plt.xlabel('Location, m')

#Moment at every point
moment_45g = forces_45g * y_tail
for i in range (1,len(y_tail)):
   moment_45g[i] = moment_45g[i-1]+moment_45g[i]
plt.subplot(143)
plt.plot(y_tail, moment_45g)
plt.ylabel('Moment, Nm')
plt.xlabel('Location, m')
plt.suptitle('Load-Shear-Moment diagrams for +9g and -4.5g static loads on tail. Location is defined from tip to root.')
plt.show()

#Deflection
deflec_45g = forces_45g * (p.b_ht/4000)**2 * (3 * p.b_ht/2000-p.b_ht/4000) /(6 * mat.E[0] * I_xx)
for i in range (2,dy+1):
    deflec_45g[dy-i] = deflec_45g[dy-i] + deflec_45g[dy-i+1]
plt.subplot(144)
plt.plot(y_tail, deflec_45g)
plt.show()
