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


########## Chord length at different spanwise locations ##########
dy = 0.001 #small spanwise section
y = np.arange(0,p.b/2., dy) #spanwise location of section
d_cLE = np.tan(p.theta_LE) * y #LE section "to be cut away from chord"
d_cTE = np.tan(p.theta_TE) * y #TE section "to be cut away from chord"

c = p.c_r - d_cLE - d_cTE #chord at each spanwise section


########## Bending at 9g ##########

CL_9g = 9. * p.g * p.MTOW / (0.5 * p.rho_0 * p.V_cruise**2 * p.S) #lift coefficient at 9g

dL_9g = CL_9g * 1./2. * p.rho_0 * p.V_cruise**2 * dy * c


moment_9g = dL_9g*y
moment_9g[0] = sum(dL_9g*y)
moment_9g_total = []

for i in range(1,len(y)):
    moment_9g[i] = moment_9g[i-1] - moment_9g[i] 
    moment_9g_total.append(moment_9g[i])

plt.plot(y[0:len(y)-1],moment_9g_total)
plt.show


