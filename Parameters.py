# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
Author: Dries
"""

# THIS FILE IS STRICTLY AND ONLY FOR PARAMETERS - NO CALCULATIONS!

import numpy as np 

Lt =  5.7 #m - Tail length
Lf = 9.04#m - fuselage length

#Mission Parameters
Nz = 9 #ultimate load factor
h_cruise = 5486.4 #5486.4m or 18000 ft - cruise altitude
V_cruise = 92.6 #m/s or 180 knots - cruise speed
rho_0 = 1.225 #density at sea level
rho_cruise = 0.698145 #density at cruise altitude
W_PL = 444 #kg - Payload
d_range = 1400 #km - Design Range
t_loiter = 2700 #s - loiter time


#Aircraft Parameters
e = 0.85 #Oswald's factor
cp = 0.000000113 #kg/J - turboprop SFC
Cd0 = 0.02 #-
n_p = 0.82 #propeller efficiency 
L_m = 0.7 #main gear length, cm
L_n = 0.7 #nose gear length, cm
MTOW = 1635 #max TO weight, in kg!!!!!!!!


#WING PARAMETERS
S = 14.31 #m2 - Wing surface area
A = 7.35 #Aspect Ratio
b = 10.35 #m - wing span
bh = 2.997 #m - horizontal wing span
c_r = 2.011 #m - root chord
c_t = 0.844 #m - tip chord
c_rh = 0.832 #m - horizontal root chord
c_th = 0.67 #m - horizontal tip chord
MAC = 1.568 #m - mean aerodynamic chord
Y = 2.327 #m - lateral MAC position
theta_LE = np.radians(3.336) 
theta_TE = np.radians(9.919)

fuel_chord_length = 0.425 #42.5% of the chord
hr1 = 0.115 #m - UPDATE IF AIRFOIL IS UPDATED
hr2 = 0.09897 #m - UPDATE IF AIRFOIL IS UPDATED
wr = fuel_chord_length * c_r #
ht1 = 0.046 #m - UPDATE IF AIRFOIL IS UPDATED
ht2 = 0.0123 #m - UPDATE IF AIRFOIL IS UPDATED
wt = fuel_chord_length * c_t
S1 = 0.5*(hr1+hr2)*wr #surface area at root
S2 = 0.5*(ht1+ht2)*wt #surface area at tip
V_fuel = Lf/3 * (S1+S2 + (S1*S2)**0.5) #maximum fuel volume in main wing

CL = 0.5

#TAIL
S_ht = 21.634 * 0.092903 #ft2 to m2
cr_ht = 0.7857
ct_ht = 0.629
b_ht = 2.828

