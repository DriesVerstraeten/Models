#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:31:38 2017

@author: driesverstraeten
"""
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

h = wm.bsection
a1 = np.abs(bs.poly1_out[2][0])+np.abs(bs.poly1_out[4][0])
c1 = np.abs(bs.poly1_out[2][-1])+np.abs(bs.poly1_out[4][-1])
a2 = c1
c2 = np.abs(bs.poly2_out[2][-1])+np.abs(bs.poly2_out[4][-1])
a3 = c2
c3 = np.abs(bs.poly3_out[2][-1])+np.abs(bs.poly3_out[4][-1])

a4 = np.abs(bs.poly1_out[3][0])+np.abs(bs.poly1_out[5][0])
c4 = np.abs(bs.poly1_out[3][-1])+np.abs(bs.poly1_out[5][-1])
a5 = c4
c5 = np.abs(bs.poly2_out[3][-1])+np.abs(bs.poly2_out[5][-1])
a6 = c5
c6 = np.abs(bs.poly3_out[3][-1])+np.abs(bs.poly3_out[5][-1])

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

print wingbox_mass()

    