#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:17:57 2017

@author: driesverstraeten
"""

import Init_Parameters as p
import Wing_model as wm
import numpy as np
import matplotlib.pyplot as plt


airfoil_coordinates = np.round(np.genfromtxt('foil1_modified.dat',skip_header=1),2)


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#PLOTTING THE AIRFOIL

plt.figure(figsize = (8.5,6),tight_layout=True)
xcoordinates = np.zeros(len(airfoil_coordinates))
ycoordinates = np.zeros(len(airfoil_coordinates))

for i in range(len(airfoil_coordinates)):
    xcoordinates[i] = airfoil_coordinates[i][0]
    ycoordinates[i] = airfoil_coordinates[i][1]

plt.plot(xcoordinates,ycoordinates)
plt.show()


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#WINGBOX SPAR LENGTHS

def wingbox_spar():
    hfront = []
    hback = []
    hfront_per = 0.15
    hback_per = 0.57
    for i in range(len(airfoil_coordinates)):
        if airfoil_coordinates[i][0] == hfront_per:
            if len(hfront) == 2:
                pass
            hfront.append(np.abs(airfoil_coordinates[i][1]))
            h_front = sum(hfront) 
        if airfoil_coordinates[i][0] == hback_per:
            if len(hback) == 2:
                pass
            hback.append(np.abs(airfoil_coordinates[i][1]))
            h_back = sum(hback)
        wingbox_width = hback_per - hfront_per
    return h_front, h_back, wingbox_width

########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#NEUTRAL AXIS LOCATION

def wingbox_NA():
    #NA = (wingbox_spar()[0]**2*p.t_skin/2 + wingbox_spar()[1]**2*p.t_skin/2 + (wingbox_spar()[3]-wingbox_spar()[2])*p.t_skin*)
    if wingbox_spar()[0] > wingbox_spar()[1]:
        h_longest = wingbox_spar()[0]
    else:
        h_longest = wingbox_spar()[1]
    return h_longest