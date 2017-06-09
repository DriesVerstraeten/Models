#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:56:43 2017

@author: driesverstraeten
"""
import MOI as mi
import Wing_model as wm
import numpy as np

f1 = mi.wingbox_MOI()[3]
f2 = mi.wingbox_MOI()[4]
y_NA = mi.wingbox_MOI()[5]
x_NA = mi.wingbox_MOI()[6]
x = mi.wingbox_MOI()[7]

def wing_spanwise_stress():
    longest = []
    for i in range(len(x)):
        longest1 = np.sqrt((f1(x)-y_NA)**2 + (x-x_NA)**2)
        longest2 = np.sqrt((f2(x)-y_NA)**2 + (x-x_NA)**2)
        longest.append(longest1)
        longest.append(longest2)
    
    print longest
    
    return longest

wing_spanwise_stress()

longest = wing_spanwise_stress()[0]
 
   
for i in range(len(longest)):
        if longest[i]==np.max(longest):
            print i