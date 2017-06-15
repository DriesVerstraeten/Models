# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 15:20:45 2017

@author: Maxime
"""
import math
import numpy as np
from Common.class12est import class2est

def calcOEW(W_fuselage,W_wing,W_tail):
#    W_avionics = 1
    OEW = W_fuselage + W_wing + W_tail # Structural weight
    
    weights,totalweight = class2est(inputlist,inputlist2,mtow,oew,fuelw,maxfuelw)
    for i in [1,2,4]:
        OEW += weights[i]
    return

#weights = np.array([W_wing,W_h,W_v,W_fuselage,W_gear,W_pwr,W_Fuel_sys,W_flightcontrol,W_iae,W_els,W_api,W_ox,W_fur])

Cd0 = 0.02
e = .85
A = 6.1
propefficiency = .82
u = .62
Vcruise = 92.6
height = 18000*0.3048
payload = 444.0
flightrange = 1400000
loitertime = 2700

aref = .614
bref = -5.887

CLmaxLand =2.1
LandingDistance = 500.

mtowinitialestimate = 1700.

sweep = 0.0
sweepquarter = math.pi/180.*sweep
taper = .3
tc = 0.15
Nult = 9.0
rootcoord = 1.76
Maxspeed = 150.0
N_pax = 4
M_D = .4
Sh = 2.24
bh = 2.97
Ah = np.sqrt(bh**2/Sh)
Ch_max = 0.833
tch = .15
trh = tch*Ch_max
lh = 6.0
Sv = 2.0
Av = 3.0
cv_max = 1.34
tcv = 0.15
Trv = tcv*cv_max
bv = np.sqrt(Av*Sv)
fuselage_l = 7.0
fuse_width = 1.2
fuse_height = 1.4
Vc = 92.6
maxfuelw = 400

sep_tanks = 2
We = 130.
Ne = 1

Nult1 = 5.7
lsm = .75
W_L = mtowinitialestimate
integralfraction = 1

inputlist = np.array([Cd0,e,A,propefficiency,u,Vcruise,height,payload,flightrange,loitertime,aref,bref,CLmaxLand,LandingDistance])
inputlist2 = np.array([sweepquarter,taper,tc,Nult,rootcoord,Maxspeed,N_pax,M_D,Sh,Ah,Ch_max,trh,bh,
                       lh,Sv,Av,Trv,bv,cv_max,fuselage_l,fuse_width,fuse_height,Vc,sep_tanks,We,Ne,
                       Nult1,lsm,W_L,integralfraction])

#mtow,oew,fuelw,WSland,weights = estweights(inputlist,inputlist2,100)