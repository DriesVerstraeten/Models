# -*- coding: utf-8 -*-
"""
Spyder Editor

Author: Dries
"""

# THIS FILE IS STRICTLY AND ONLY FOR PARAMETERS - NO CALCULATIONS!

import numpy as np 


#ATMOSPHERIC PARAMETERS
rho_0 = 1.225 #density at sea level
rho_cruise = 0.698145 #density at cruise altitude
g = 9.80665 #gravitational acceleration


#MISSION PARAMETERS
Nz = 12. #ultimate load factor
h_cruise = 7620 #5486.4m or 18000 ft - cruise altitude
V_cruise = 123.5 #m/s or 180 knots - cruise speed
W_PL = 444. #kg - Payload
d_range = 1400. #km - Design Range
t_loiter = 2700. #s - loiter time


#AIRCRAFT PARAMETERS
L_m = 0.7 #main gear length, cm
L_n = 0.7 #nose gear length, cm
MTOW = 1712.16 #max TO weight, in kg!!!!!!!!

#FUSELAGE PARAMETERS
Lt =  4.5 #m - Tail length
Lf = 10   #m - fuselage length


#WING PARAMETERS
S = 15.4465 #m2 - Wing surface area
A = 7.35 #Aspect Ratio
b = 10.655 #m - wing span
#bh = 2.997 #m - horizontal wing span
#c_r = 2.011 #m - root chord
#c_t = 0.844 #m - tip chord
#c_rh = 0.832 #m - horizontal root chord
#c_th = 0.67 #m - horizontal tip chord
MAC = 1.5746 #m - mean aerodynamic chord
#Y = 2.327 #m - lateral MAC position
X_le_mac = 2.04
theta_LE = np.radians(3.336) 
theta_TE = np.radians(9.919)

hr1 = 0.115 #m - UPDATE IF AIRFOIL IS UPDATED
hr2 = 0.09897 #m - UPDATE IF AIRFOIL IS UPDATED
ht1 = 0.046 #m - UPDATE IF AIRFOIL IS UPDATED
ht2 = 0.0123 #m - UPDATE IF AIRFOIL IS UPDATED
#wt = fuel_chord_length * c_t
#S1 = 0.5*(hr1+hr2)*wr #surface area at root
#S2 = 0.5*(ht1+ht2)*wt #surface area at tip
#V_fuel = Lf/3. * (S1+S2 + (S1*S2)**0.5) #maximum fuel volume in main wing

t_skin = 0.002 #m - skin thickness


#TAIL
S_ht = 1.7756 #ft2 to m2
cr_ht = 1.44959 #root cord hor tail
ct_ht = 0.3*cr_ht #tip cord hor tail
b_ht = 1.88446 #span hor tail
#W_ht = 13.432 #mass hor tail


#MATERIAL PROPERTIES
rho_CF = 1800 #kg/m3
E_CF = 365*10**9 #Pa, 365GPa
sigma_tensile = 4500*10**6 #MPa


