#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:08:10 2017

@author: driesverstraeten
"""

#Calculations for the WING MODEL

import Parameters as p
import numpy as np
import matplotlib.pyplot as plt
#lalal
plt.close()
fig = plt.figure(figsize = (8.5,6),tight_layout=True)

########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#CHORD LENGTH AT DIFFERENT SPANWISE LOCATIONS

dy = 0.001 #small spanwise section
y = np.arange(0,p.b/2.+dy, dy) #spanwise location of section
y1 = np.arange(0,p.b/2.+2*dy, dy)
d_cLE = np.tan(p.theta_LE) * y #LE section "to be cut away from chord"
d_cTE = np.tan(p.theta_TE) * y #TE section "to be cut away from chord"

c = p.c_r - d_cLE - d_cTE #chord at each spanwise section

CL_9g = 9. * p.g * p.MTOW / (0.5 * p.rho_0 * p.V_cruise**2. * p.S) #lift coefficient at 9g
CL_45g = -4.5 * p.g * p.MTOW / (0.5 * p.rho_0 * p.V_cruise**2. * p.S) #Lift coefficient at -4.5g
                            
########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#SHEAR AT 9G

dL_9g = CL_9g * 1./2. * p.rho_0 * p.V_cruise**2. * dy * c #the small lift contribution from every section


dL_9g_total = np.zeros(len(y)+1) #make a list of zeroes to later overwrite these in the next loop
          
for i in range(0,len(y)+1):
    if i == 0:
        dL_9g_total[i] = sum(CL_9g * 1./2. * p.rho_0 * p.V_cruise**2. * dy * c) #overwrite the zeroes
    elif i == 1:
        dL_9g_total[i] = dL_9g_total[i-1] - dL_9g[0]
    else:
        dL_9g_total[i] = dL_9g_total[i-1] - dL_9g[i-1]

ax1 = fig.add_subplot(221)
ax1.plot(y1,dL_9g_total)  
ax1.set_title('Shear force at 9g')
ax1.set_ylabel('Shear force [N]')
ax1.set_xlabel('Wing span [m]')
ax1.set_ylim(dL_9g_total[-1],dL_9g_total[0])
ax1.set_xlim([y[0],y[-1]])


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#BENDING AT 9g

dM_9g = dL_9g * y

dM_9g_total = np.zeros(len(y)+1)

for i in range(0,len(y)+1):
    if i == 0:
        dM_9g_total[i] = sum(dM_9g)
    elif i == 1:
        dM_9g_total[i] = dM_9g_total[i-1] - dM_9g[0]
    else:
        dM_9g_total[i] = dM_9g_total[i-1] - dM_9g[i-1]


ax2 = fig.add_subplot(222)
ax2.plot(y1,dM_9g_total)  
ax2.set_title('Bending moment at 9g')
ax2.set_ylabel('Bending moment [Nm]')
ax2.set_xlabel('Wing span [m]')
ax2.set_ylim(dM_9g_total[-1],dM_9g_total[0])
ax2.set_xlim([y[0],y[-1]])

########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#BENDING AT -4.5g

dL_45g = CL_45g * 1./2. * p.rho_0 * p.V_cruise**2. * dy * c #the small lift contribution from every section


dL_45g_total = np.zeros(len(y)+1) #make a list of zeroes to later overwrite these in the next loop
          
for i in range(0,len(y)+1):
    if i == 0:
        dL_45g_total[i] = sum(CL_45g * 1./2. * p.rho_0 * p.V_cruise**2. * dy * c) #overwrite the zeroes
    elif i == 1:
        dL_45g_total[i] = dL_45g_total[i-1] - dL_45g[0]
    else:
        dL_45g_total[i] = dL_45g_total[i-1] - dL_45g[i-1]

ax3 = fig.add_subplot(223)
ax3.plot(y1,dL_45g_total)  
ax3.set_title('Shear force at -4.5g')
ax3.set_ylabel('Shear force [N]')
ax3.set_xlabel('Wing span [m]')
ax3.set_ylim(dL_45g_total[-1],dL_45g_total[0])
ax3.set_xlim([y[0],y[-1]])


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#BENDING AT -4.5g

dM_45g = dL_45g * y

dM_45g_total = np.zeros(len(y)+1)

for i in range(0,len(y)+1):
    if i == 0:
        dM_45g_total[i] = sum(dM_45g) #the first point is the total moment
    elif i == 1:
        dM_45g_total[i] = dM_45g_total[i-1] - dM_45g[0] #this is the second point, this is so that the last point is equal to 0
    else:
        dM_45g_total[i] = dM_45g_total[i-1] - dM_45g[i-1] #all the other points

ax4 = fig.add_subplot(224)
ax4.plot(y1,dM_45g_total)  
ax4.set_title('Bending moment at -4.5g')
ax4.set_ylabel('Bending moment [Nm]')
ax4.set_xlabel('Wing span [m]')
ax4.set_ylim([dM_45g_total[0],dM_45g_total[-1]])
ax4.set_xlim([y[0],y[-1]])











plt.show()