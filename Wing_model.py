#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:08:10 2017

@author: driesverstraeten
"""

#Calculations for the WING MODEL

import Parameters as p
import numpy as np


########## Chord length at different spanwise locations ##########
dy = 0.01 #small spanwise section
y = np.arange(0,p.b/2., dy) #spanwise location of section
d_cLE = np.tan(p.theta_LE) * y #LE section "to be cut away from chord"
d_cTE = np.tan(p.theta_TE) * y #TE section "to be cut away from chord"

c = p.c_r - d_cLE - d_cTE #chord at each spanwise section


########## Shear at 9g ##########

CL_9g = 9. * p.g * p.MTOW / (0.5 * p.rho_0 * p.V_cruise**2 * p.S) #lift coefficient at 9g

dL_9g = CL_9g * 1./2. * p.rho_0 * p.V_cruise**2 * dy * c


moment = dL_9g*y
moment[0] = sum(dL_9g*y)


for i in range(1,len(y)):
    moment[i] = moment[i-1] - moment[i] 
    print moment[i]




dM_9g = dL_9g * y

