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
dy = 0.01
y = np.arange(0,p.b/2, dy) #spanwise section
d_cLE = np.tan(p.theta_LE) * y #LE section "to be cut away from chord"
d_cTE = np.tan(p.theta_TE) * y #TE section "to be cut away from chord"

c = p.c_r - d_cLE - d_cTE #chord at each spanwise section


########## Bending Diagram ##########

dL = 1./2. * p.rho_0 * p.V_cruise**2 * dy * c * p.CL
L = np.sum(dL)
