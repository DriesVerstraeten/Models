#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 09:30:19 2017

@author: driesverstraeten
"""

import Wing_model as wm
import Init_Parameters as p

c_r = wm.C1
c_t = wm.C4
S = p.S
h = p.b/2

c_hr = p.cr_ht
c_ht = p.ct_ht
Sh = p.S_ht
ht = 1.71

L_tail = p.Lt

Iyy_ht = ht * (c_ht + c_hr) * (c_ht**2 + 7 * c_hr**2)/48 * 2 + Sh * L_tail **2
Iyy_w = h * (c_t + c_r) * (c_t**2 + 7 * c_r**2)/48 * 2 
Iyy_tot = Iyy_ht + Iyy_w

Ixx_ht = ht**3*(3*c_ht+c_hr)/12 * 2 
Ixx_w = h**3*(3*c_t+c_r)/12 * 2 
Ixx_tot = Ixx_ht + Ixx_w
               
Izz_ht = p.MAC*(p.MAC*0.14)**3/12 - (p.MAC-p.t)*((p.MAC-p.t)*0.14)**3/12
Izz_w = p.MAC*(p.MAC*0.14)**3/12 - (p.MAC-p.t)*((p.MAC-p.t)*0.14)**3/12
 