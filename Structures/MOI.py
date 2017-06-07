#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:17:57 2017

@author: driesverstraeten
"""

import Init_Parameters as p
import Wing_model as wm
import numpy as np

airfoil_coordinates = np.round(np.genfromtxt('NACA63215.txt',skip_header=1),2)

def wingbox_spar():
    h15 = []
    h55 = []
    for i in range(len(airfoil_coordinates)):
        if airfoil_coordinates[i][0] == 0.15:
            h15.append(np.abs(airfoil_coordinates[i][1]))
            h_15 = sum(h15)  
    for i in range(len(airfoil_coordinates)):
        if airfoil_coordinates[i][0] == 0.55:
            #print airfoil_coordinates[i][1]
            h55.append(np.abs(airfoil_coordinates[i][1]))
            h_55 = sum(h55)
    
    
    return 

