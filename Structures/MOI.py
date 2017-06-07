#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:17:57 2017

@author: driesverstraeten
"""

import Init_Parameters as p
import Wing_model as wm
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt


t = p.t_skin

airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#WINGBOX SPAR LENGTHS

def wingbox_spar():
    
    hfront = []
    hback = []
    hfront_per = 0.15
    hback_per = 0.57
    y_front_bottom = []
    y_back_bottom = []
    
    airfoil_coordinates = np.round(np.genfromtxt('foil1_modified.dat',skip_header=1),2)
    
    for i in range(len(airfoil_coordinates)):
        if airfoil_coordinates[i][0] == hfront_per:
            if len(hfront) == 2:
                continue
            hfront.append(np.abs(airfoil_coordinates[i][1]))
            h_front = sum(hfront)
            y_front_bottom.append(airfoil_coordinates[i][1])
            
        if airfoil_coordinates[i][0] == hback_per:
            if len(hback) == 2:
                continue
            hback.append(np.abs(airfoil_coordinates[i][1]))
            h_back = sum(hback)
            y_back_bottom.append(airfoil_coordinates[i][1])

        wingbox_width = hback_per - hfront_per
        #difference_spar = np.abs(h_back-h_front)
        
    return h_front, h_back, y_front_bottom, y_back_bottom, wingbox_width,


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#NEUTRAL AXIS LOCATION

def wingbox_NA():
    
    if wingbox_spar()[0] > wingbox_spar()[1]:
        h_longest = wingbox_spar()[0]
        h_shortest = wingbox_spar()[1]
    else:
        h_longest = wingbox_spar()[1]
        h_shortest = wingbox_spar()[0]
        
    y_datum = -0.07
    
    front_spar_area = wingbox_spar()[0]*t
    back_spar_area = wingbox_spar()[1]*t
    
    front_spar_cg = wingbox_spar()[0]/2 + np.abs(y_datum - wingbox_spar()[2][1])
    back_spar_cg = wingbox_spar()[1]/2 + np.abs(y_datum - wingbox_spar()[3][1])
    
    bottom_skin_width = 0
    
    print back_spar_cg
    
    #NA = (wingbox_spar()[0]**2*p.t_skin/2 + wingbox_spar()[1]**2*p.t_skin/2 + wingbox_spar()[2]*p.tskin*(h_shortest+wingbox_spar[3]/2)/(wingbox_spar()[0]*p.t_skin + wingbox_spar()[1]*p.t_skin + wingbox_spar()[2]*p.tskin*(h_shortest+wingbox_spar[3]/2
    
    return h_longest

wingbox_NA()


"""
########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#PLOTTING THE AIRFOIL

plt.figure(figsize = (8.5,6),tight_layout=True)
xcoordinates = np.zeros(len(airfoil_coordinates))
ycoordinates = np.zeros(len(airfoil_coordinates))

for i in range(len(airfoil_coordinates)):
    xcoordinates[i] = airfoil_coordinates[i][0]
    ycoordinates[i] = airfoil_coordinates[i][1]

#f = interpolate.interp1d(xcoordinates, ycoordinates, kind='linear')
#print xcoordinates, ycoordinates

xfront = 0.15 * np.ones(2)
xback = 0.57 * np.ones(2)
yfront = [-0.04,0.07]
yback = [-0.04,0.09]

plt.plot(xfront,yfront)
plt.plot(xback,yback)


plt.plot(xcoordinates,ycoordinates)


plt.show()"""























