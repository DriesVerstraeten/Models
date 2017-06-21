#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:31:38 2017

@author: driesverstraeten
"""

import time
start_time = time.time()

import MOI as mi
import Init_Parameters as p
import Wing_model as wm
import wing_bending_stress as bs
import numpy as np



FS_location = mi.FS_location
BS_location = mi.BS_location
c = wm.c
dy = wm.dy
t = p.t_skin

poly1_out = bs.poly1_out
poly2_out = bs.poly2_out
poly3_out = bs.poly3_out

h = wm.bsection
a1 = np.abs(poly1_out[2][0])+np.abs(poly1_out[4][0])
c1 = np.abs(poly1_out[2][-1])+np.abs(poly1_out[4][-1])
a2 = c1
c2 = np.abs(poly2_out[2][-1])+np.abs(poly2_out[4][-1])
a3 = c2
c3 = np.abs(poly3_out[2][-1])+np.abs(poly3_out[4][-1])

a4 = np.abs(poly1_out[3][0])+np.abs(poly1_out[5][0])
c4 = np.abs(poly1_out[3][-1])+np.abs(poly1_out[5][-1])
a5 = c4
c5 = np.abs(poly2_out[3][-1])+np.abs(poly2_out[5][-1])
a6 = c5
c6 = np.abs(poly3_out[3][-1])+np.abs(poly3_out[5][-1])

a_FS = np.hstack((a1,a2,a3))
c_FS = np.hstack((c1,c2,c3))

def wingbox_area(): #this is for the entire wing, not just half like in the other files!
    area_wingbox_skin = (BS_location - FS_location)*sum(c*dy)*2*2 #area of the upper and lower skin of wingbox
    area_FS = sum(h*(a_FS+c_FS)/2) * 2 #multiplied by 2 for the entire wing
    area_BS = sum(h*(a_FS+c_FS)/2) * 2
                        
    return area_wingbox_skin, area_FS, area_BS

area_wingbox_skin, area_FS, area_BS = wingbox_area()

def wingbox_volume():
    volume_wingbox_skin = area_wingbox_skin * 2 * t 
    volume_FS = area_FS * t 
    volume_BS = area_BS * t
    total_volume = volume_wingbox_skin + volume_FS + volume_BS
    
    return total_volume, volume_wingbox_skin, volume_FS, volume_BS

total_volume, volume_wingbox_skin, volume_FS, volume_BS = wingbox_volume()

def wingbox_mass():
    wingbox_mass = total_volume * p.rho_CF
    return wingbox_mass

wingbox_mass = wingbox_mass()

Ixx, Iyy, Ixy, x_NA, y_NA, x_span, LE_area, TE_area, ribs_area, total_volume_wings = mi.wingbox_MOI_total()

def LE_TE_mass():
    LE_volume = 2* LE_area * t/2 
    TE_volume = 2* TE_area * t /2
    volume_total = LE_volume + TE_volume
    mass = volume_total * p.rho_CF
    return mass

LE_TE_mass = LE_TE_mass()

def rib_mass():
    ribs_volume = 2* ribs_area * t
    ribs_mass = sum(ribs_volume * p.rho_CF)
    return ribs_mass
    
rib_mass = rib_mass()
    
def wing_total_mass():
    return wingbox_mass + LE_TE_mass + rib_mass

print "WEIGHT CALC:", wing_total_mass(), "kg"
    
    
    
print("WEIGHT CALCULATOR --- %s seconds ---" % (time.time() - start_time))

    