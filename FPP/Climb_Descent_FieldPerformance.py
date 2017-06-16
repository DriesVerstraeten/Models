#NOTES
#
#Power should still be implemented as a function of altitude
#Power used is mainly in BHP
#P_TO = P_BHP at sea level
#
#


#Importing modules
import numpy as np
import matplotlib.pyplot as plt
import os

#Importing functions
import Common.CalcISA as ISA
from FPP.propulsion import Analyse_prop

#ISA lists
Ttab = []
Rhotab = []
htab = []
ptab = []
State = 3*[0]
for h in range(0,84000):
    State = ISA.ISA(h)
    Ttab.append(State[0])
    Rhotab.append(State[2])
    ptab.append(State[1])
    htab.append(h)


airfoil_path = os.getcwd() + '\Polars\Eppler_prop.npz'
rpm = 2800

###############################################################################
##############################  INPUTS  #######################################
###############################################################################

#Initial test input values (IMPERIAL UNITS)

#CL_0 = 0.39
#CL_alpha = 0.1
#alpha_TO = 2. #deg
#MTOW = 3400. #lbf
#rho = 0.002378 #slugs/ft3
#S = 144.9 #ft2
#CL_max_TO = 2.2 
#etha_p = 0.85 #propeller efficiency
#P_TO = 450. #BHP
#V_c = 287. #ft/s - Cruising velocity
#V_H = 314. #ft/s - Maximum velocity
#D_prop = 5.577427822 #ft
#D_spinner = 1.312335958 #ft

#CD_0 = 0.0250 
#A = 7 #Aspect Ratio
#e = 0.85 #Oswald efficiency factor
#LD_max = 16. #Maximum lift-over-drag


###############################################################################
#########################  Take-off performance  ##############################
###############################################################################

"""
## Inputs ##

# MTOW [N] Take-off weight, usually MTOW, but can be adapted for aerobatic weight
# S [m2] Surface area
# A [-] Aspect ratio
# e [-] Oswald efficiency factor
# CD0 [-] Zero lift drag coefficient of the wing
# etha_p [-] Prop efficiency
# P_TO [W] Take-off power
# alpha_TO [deg] angle of attack during ground run
# CL_max_TO [-] Maximum lift coefficent in take-off configuration
# V_c [m/s] Cruise speed
# V_H [m/s] Maximum speed


# CAN BE REPLACED BY CL_TO_rw
# CL_0 [-] Lift coefficient at zero angle of attack
# CL_alpha [-] dCL/dalpha

# CAN BE REPLACED BY Tstatic
# D_prop [-] Propellor diameter
# D_spinner [-] Spinner diameter


## Outputs ##

# CL_TO [-] Lift coefficient during ground run part of take-off
# CD_TO [-] Drag coefficient during ground run part of take-off
# V_S1_si [m/s] Stall speed in landing config 
# V_LOF_si [m/s] Lift off speed 
# S_TO_si [m] Take-off distance

#Verification test PC-12 
#Takeoff(46500, 25.81, 10.27, 0.85, 0.04, 0.8, 895, 0.4, 0.1, 2, 2.1, 2.67, 0.3, 146, 155)
#Take-off distance = 810m (model) vs 793m (reality)
#It should be noted that the CL, CD, prop efficiency and e are a rough estimation however

#Cirrus Takeoff(15123.95, 13.387, 10, 0.757, 0.0350, 0.85, 231.167, 0.39, 0.1, 2, 1.69, 1.93, 0.4572, 87.46, 95.7)
#Gives 488 m (model) vs. 468 m (book)
#Difference is explained by higher vstall than calculated
"""

def Takeoff(MTOW, S, A, e, CD_0, CL_0, CL_alpha, alpha_TO, CL_max_TO, V_c, V_H):


    #NOTES:
    #All input values and formulas are performed using imperial units
    #The method used is described in chapter 17 method #3
    
    #Importing power parameters
    T_static, NA1, NA2, P_TO, NA3 = Analyse_prop(airfoil_path, 0, 0, rpm)
    
    #Weight from kg to N
    MTOW = MTOW*9.80665
    
    # Not an input for the function, however an input locally
    rho_sealevel = 0.002378 #slugs/ft3
    mu_g_nobrakes = 0.08 #Ground friction coefficient
    g = 32.1740 #ft/s2
    h_obst = 50. #ft
    
    #Converting the input values from si units to imperial units
    MTOW = MTOW*0.224808943871 #Newton to lbf
    S = S*10.7639104 #m2 to ft2
    #D_prop = D_prop*3.2808399 #m to ft
    #D_spinner = D_spinner*3.2808399 #m to ft
    V_c = V_c*1.943844 #m/s to kts
    V_H = V_H*1.943844 #m/s to kts
    P_TO = P_TO*1.3410220888/1000 #W to BHP
    
    #Basic functions taken from Chapter 17 of the General Aviation book
    CL_TO = CL_0 + CL_alpha * alpha_TO 
    CD_TO = CD_0 + ((CL_TO**2)/(np.pi*A*e))
    V_S1 = np.sqrt((2.*MTOW)/(rho_sealevel*S*CL_max_TO)) #ft/s
    V_LOF = 1.1 * V_S1 #ft/s
    
    V_LOF_ms = 0.3048*V_LOF
    BLA1, BLA2, BLA3, BLA4, etha_p = Analyse_prop(airfoil_path, 0, V_LOF_ms, rpm)
    
    T_c = (etha_p*550.*P_TO)/(V_c*1.68781) #lbf
    T_H = (etha_p*550.*P_TO)/(V_H*1.68781) #lbf  
    
  
    #Setting up the factors of a cubic spline as is explained on page 810 from the book
    Q = np.matrix([[0., 0., 0., 1.],
                   [(V_c**3.), (V_c**2.), V_c, 1.],
                   [(3.*V_c**2.), (2.*V_c), 1., 0.],
                   [(V_H**3.), (V_H**2.), V_H, 1.]])
    W = np.matrix([[T_static],
                   [T_c],
                   [(-etha_p*325.8*P_TO)/(V_c**2.)],
                   [T_H]])
    x = np.linalg.solve(Q, W)
    O = float(x[0])
    B = float(x[1])
    C = float(x[2])
    E = float(x[3])
    
    
    i = 0
    time = 0
    dt = 1
    t = []
    V = [0.]
    L = []
    D = []
    q = []
    a = []
    S1 = []
    S2 = []
    S_G = []
    T_V = []
    # Pas terug aan naar 10000
    
    for i in range(0, 1):
        if i == 0:
            time = 0
            thrust = E
            dynam_pressure = 0.5*rho_sealevel*V[i]**2
            lift = dynam_pressure*S*CL_TO
            drag = dynam_pressure*S*CD_TO
            acceleration = (thrust-drag-(mu_g_nobrakes*(MTOW - lift)))/(MTOW/g)
            distance_1 = V[i] * dt
            distance_2 = 0.
            distance_total = distance_1 + distance_2
            
            t.append(time)        
            T_V.append(thrust)
            q.append(dynam_pressure)
            L.append(lift)
            D.append(drag)
            a.append(acceleration)
            S1.append(distance_1)
            S2.append(distance_2)
            S_G.append(distance_total)
        
        if i==1:
            time = time + dt
            thrust = (O*(V[0]*0.592483801)**3.) + (B*(V[0]*0.592483801)**2.) + (C*(V[0]*0.592483801)) + E
            dynam_pressure = 0.5*rho_sealevel*V[i-1]**2
            lift = dynam_pressure*S*CL_TO
            drag = dynam_pressure*S*CD_TO
            acceleration = (thrust-drag-(mu_g_nobrakes*(MTOW - lift)))/(MTOW/g)
            airspeed = V[i-1] + acceleration*dt
            distance_1 = V[0] * dt
            distance_2 = 0.5 * a[0] * (dt**2)
            distance_total = S_G[0] + distance_1 + distance_2
            
            t.append(time)        
            T_V.append(thrust)
            q.append(dynam_pressure)
            L.append(lift)
            D.append(drag)
            a.append(acceleration)
            V.append(airspeed)
            S1.append(distance_1)
            S2.append(distance_2)
            S_G.append(distance_total)
        
        if i>1:
            time = time + dt
            thrust = (O*(V[i-1]*0.592483801)**3.) + (B*(V[i-1]*0.592483801)**2.) + (C*V[i-1]*0.592483801) + E
            dynam_pressure = 0.5*rho_sealevel*V[i-1]**2
            lift = dynam_pressure*S*CL_TO
            drag = dynam_pressure*S*CD_TO
            acceleration = (thrust-drag-(mu_g_nobrakes*(MTOW - lift)))/(MTOW/g)
            airspeed = V[i-1] + acceleration*dt
            distance_1 = V[i-1] * dt
            distance_2 = 0.5 * a[i-1] * (dt**2)
            distance_total = S_G[i-1] + distance_1 + distance_2
            
            t.append(time)        
            T_V.append(thrust)
            q.append(dynam_pressure)
            L.append(lift)
            D.append(drag)
            a.append(acceleration)
            V.append(airspeed)
            S1.append(distance_1)
            S2.append(distance_2)
            S_G.append(distance_total)
            
            
    S_ROT = V_LOF
    VLOF_list = [V_LOF]*len(S_G) #Only useful for plotting
    
    V_TR = 1.15*V_S1*0.592483801 #kts
    V_TR_ft = 1.15*V_S1 #ft/s
    
    loc = np.where(np.array(V) < V_TR_ft) #V in kts
    loc_list = loc[0][-1]
    
    loc_VLOF = np.where(np.array(V) < V_LOF) #V in ft/s
    loc_list_VLOF = loc_VLOF[0][-1]
    
    CL_VTR = (2*MTOW)/(rho_sealevel*(V_TR_ft**2)*S) #V in ft/s
    k = 0.04207#1./(pi*A*e)
    CD_VTR = CD_0 + k*(CL_VTR**2)
    
    S_TR = (0.2156*((V_S1)**2)*((T_V[loc_list]/MTOW)-(1/(CL_VTR/CD_VTR)))) #ft
    
    theta_climb = np.arcsin(((T_V[loc_list]/MTOW)-(1/(CL_VTR/CD_VTR))))*(180/np.pi) #degrees
    R_transition = (0.2156*(V_S1)**2) #ft
    h_transition = R_transition*(1. - np.cos(theta_climb*(np.pi/180))) #ft
    
    h_obst = 50 #ft, this is a CS23 requirement
    S_C = (h_obst-h_transition)/(np.tan(theta_climb*(np.pi/180))) #ft
    
    #SI units outputs
    S_TO_si = (S_G[loc_list_VLOF] + S_ROT + S_TR + S_C)*0.3048
    V_LOF_si = V_LOF*0.3048
    V_S1_si = V_S1*0.3048

      
    return(CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si, T_static, P_TO, etha_p)


###############################################################################
########################  Climb performance  ##################################
###############################################################################

"""
## Inputs ##

# MTOW [N] Take-off weight, usually MTOW, but can be adapted for aerobatic weight
# S [m2] Surface area
# A [-] Aspect ratio
# e [-] Oswald efficiency factor
# CD0 [-] Zero lift drag coefficient of the wing
# etha_p [-] Prop efficiency
# P [W] Power during climb --------------------------------> Dependent on altitude
# LD_max [-] Maximum lift-over-drag ratio
# hcruise [m] Cruise altitude

## Outputs ##

# t_cruiseh [s] Time it takes to get to cruise altitude
# ROC_SL_si [m/s] Maximum Rate of Climb at sea level
# h_absolute_si [m] Absolute service ceiling, meaning the altitude at which the ROC = 0

# Verification using PC-12
# Climb(46500, 25.81, 10.27, 0.85, 0.04, 0.8, 800, 12, 9150)
# Climb rate = 9.9 m/s (model) vs. 9.75 m/s (flight test)
"""
#Minimum drag coefficient instead of Cd0
def Climb(MTOW, S, A, e, CD_0, etha_p, P, LD_max, hcruise):
    
    #Weight from kg to N
    MTOW = MTOW*9.80665
    
    #Converting the input values from si units to imperial units
    MTOW = MTOW*0.224808943871 #Newton to lbf    
    S = S*10.7639104 #m2 to ft2
    P = P*1.3410220888*1000 #kW to BHP
    

    #Airspeed for best ROC
    V_Y_list = []
    hmax = 28000*0.3048 #meter
    h0 = 0.
   
    
    for i in range(0, hmax):
        rho_h = Rhotab[i] * 0.00194032033 #Converted to slugs/ft3       
        V_Y = np.sqrt((2./rho_h)*(MTOW/S)*np.sqrt((1./(np.pi*A*e))/(3.*CD_0)))
        V_Y_list.append(V_Y)
    
    #Rate of climb
    ROC_list = []
    
    #If power changes according to height try changing the V_Y to a constant value
    
    for i in range(0, len(V_Y_list)):
        ROC = 60.*(((etha_p*550.*P)/MTOW)-(V_Y_list[i]*(1.1547/LD_max))) #ft/min, CHANGE P_TO to power dependent on altitude
        ROC_si = ROC*0.00508
        ROC_list.append(ROC_si)   
        ROC_SL_si = ROC_list[0]
        #ROC_SL_si = ROC_SL*0.00508
        if ROC<0.2 and ROC>-0.2:
            h_absolute_si = htab[i]
        if ROC<100.1 and ROC>99.9:
            h_service_si = htab[i]
        
    #It should be noted that the maximum climb rate is calculated at MTOW with full fuel tanks
    
    #Polynomial fit for the ROC vs. Altitude plot
    z = np.polyfit(htab[0:len(ROC_list)], ROC_list, 3) #Third order
    polyfit = []
    for i in range(0, len(V_Y_list)):
        Z = z[0]*(htab[i]**3) + z[1]*htab[i]**2 + z[2]*htab[i] + z[3] #Third order
        #Z = z[0]*htab[i] + z[1] #First order
        polyfit.append(Z)
    
    
    #polynomial fit
    aa = z[0]
    bb = z[1]
    cc = z[2]
    dd = z[3]
    
    
    #Time to altitiude & service ceiling
    #ROC_a = (aa*(hcruise-h0))/(np.log(aa*hcruise + bb) - np.log(aa*h0 + bb))
    #t = (hcruise - h0)/ROC_a
    
    #Modified formula from book page 837
    ROC_a_3rd = (aa*((hcruise-h0)**3)+(bb*((hcruise-h0)**2)+(cc*(hcruise-h0))))/(np.log((aa*(hcruise**3))+(bb*(hcruise**2))+(cc*hcruise) + dd) - np.log((aa*(h0**3))+(bb*(h0**2))+(cc*h0) + dd))
    t_cruiseh = (hcruise - h0)/ROC_a_3rd
    
    return(t_cruiseh, ROC_SL_si, h_absolute_si)


###############################################################################
######################  Descent performance  ##################################
###############################################################################

"""
## Inputs ##

# W_D [N] The weight of the aircraft right after it reached cruise altitude (not MTOW)
# S [m2] Surface area
# A [-] Aspect ratio
# e [-] Oswald efficiency factor
# CD_0 [-] Zero lift drag coefficient
# LD_max [-] Maximum lift-over-drag ratio
# hcruise [m] Cruise altitude


## Outputs ##

# R_glide [m] The glide distance starting from cruise altitude and under optimal glide conditions

# Verification using PC-12 flight test data
# Descent(46500, 25.81, 10.27, 0.85, 0.04, 12, 9150)
# gives 110 km, flight test gave 140 km
# https://www.flightglobal.com/news/articles/flight-test-pilatus-pc-12-power-of-one-187732/
# https://en.wikipedia.org/wiki/Pilatus_PC-12
"""

def Descent(MTOW, S, A, e, CD_0, LD_max, hcruise):
    
    #Weight from kg to N
    MTOW = MTOW*9.80665
    
    W_descent = 0.965*MTOW
    
    #Best glide speed
    V_BG_list = [] #ft/s
    
    for i in range(0, hcruise):
        rho_h = Rhotab[i] * 0.00194032033 #Converted to slugs/ft3       
        V_BG = np.sqrt((2/rho_h)*(np.sqrt((1./(np.pi*A*e))/CD_0))*(W_descent/S)) #What weight are we going to take?
        V_BG_list.append(V_BG)
    
    V_BG_list = V_BG_list[::-1] #Reversing the list, we're descending here
    
    #Minimum glide angle
    theta_min_glide = np.arctan(1./LD_max)
    
    
    #Glide distance - provided the best glide speed and minimum glide angle
    R_glide = hcruise * LD_max
    
    return(R_glide) 

###############################################################################
######################  Landing performance  ##################################
###############################################################################
"""
## inputs ##

# CL_landing_max [-] Maximum lift coefficient in landing configuration
# T_static [N] Static thrust of the engine --------------------------------> Or just power input
# theta_app [deg] approach angle --------------------------> CS23??
# CL_landing_td [-] Lift coefficient after toachdown -------------> Same as during ground run on TO?
# A [-] Aspect ratio
# e [-] Oswald efficiency factor
# etha_p_landing [-] Propellor efficiency during landing
# W_landing [N] Landing weight ------------------------------> What this is should be specified
# S [m2] Surface area


## outputs ##

# S_landing_si [m] Total landing distance
# S_ground_roll_si [m] Ground roll distance
# S_landing_reverse_si [m] Total landing distance using reverse thrust
# S_ground_roll_reverse_si [m] Ground roll distance using reverse thrust
# V_REF [m/s] Approach landing speed

# Verification using PC-12
# Landing(44145, 25.81, 10.27, 0.85, 0.04, 0.8, 2.1, 8000, 3, 0.6)
# Assumed was a Tstatic of 8000 N, CL_max_landing of 2.1 and a CL after touchdown of 0.6
# Landing distance = 677 m (Model) vs. 661 m (Flight test) on dry ashphalt (mu_g = 0.3)
"""

def Landing(MTOW, S, A, e, CD_0, etha_p_landing, CL_landing_max, T_static, CL_landing_td):
    
    #Weight from kg to N
    MTOW = MTOW*9.80665
    
    #Defined the weight
    W_landing = MTOW
    
    # Not an input for the function, however an input locally
    
    rho_sealevel = 0.002378 #slugs/ft3
    mu_g_brakes = 0.2 #Ground friction coefficient
    g = 32.1740 #ft/s2
    h_obst = 50 #ft
    
    #Converting 
    
    S = S*10.7639104 #m2 to ft2
    T_static = T_static*0.224808943871 #Newton to lbf
    W_landing = W_landing*0.224808943871 #Newton to lbf

    # Functions
    
    CD_landing_td = CD_0 + ((CL_landing_td**2)/(np.pi*A*e))
    
    T_landing = 0.07*T_static #lbf
    T_landing_reverse_thrust = -0.6*T_static #lbf

    V_SO = np.sqrt((W_landing/S)*(2./rho_sealevel)*(1./CL_landing_max)) #Landing stall speed in ft/s
    V_REF = 1.3*V_SO #Approach speed
    V_FLR = 1.3*V_SO #Flare speed in ft/s
    V_TD = 1.1 * V_SO #Touchdown speed in ft/s
    V_BR = 1.1 * V_SO #Brake speed in ft/s
    
    P_BHP_idle = T_landing*V_TD * 0.001818181818179 #Converts lbf*ft/s to horsepower
    P_BHP_reverse = T_landing_reverse_thrust*V_TD * 0.001818181818179 #Converts lbf*ft/s to horsepower

    L_landing_td = 0.5*rho_sealevel*((V_BR/np.sqrt(2))**2)*S*CL_landing_td
    D_landing_td = 0.5*rho_sealevel*((V_BR/np.sqrt(2))**2)*S*CD_landing_td
    
    theta_app = np.arcsin((1./(CL_landing_td/CD_landing_td))-(T_landing/W_landing))*(180./np.pi)

    h_f = 0.1512*(V_SO**2.)*(1.-np.cos(theta_app*(np.pi/180.)))
    S_A = (h_obst - h_f)/(np.tan(theta_app*(np.pi/180.)))
    S_F = 0.1512*(V_SO**2.)*np.sin(theta_app*(np.pi/180))
    S_FR = V_TD
    S_BR = abs(((V_BR**2.)*W_landing)/(2.*g*((np.sqrt(2.)*etha_p_landing*550.*(P_BHP_idle/V_BR))-D_landing_td-(mu_g_brakes*(W_landing-L_landing_td)))))
    S_BR_reverse = abs(((V_BR**2.)*W_landing)/(2.*g*((np.sqrt(2.)*etha_p_landing*550.*(P_BHP_reverse/V_BR))-D_landing_td-(mu_g_brakes*(W_landing-L_landing_td)))))

    S_landing  = S_A + S_F + S_FR + S_BR
    S_ground_roll = S_FR + S_BR
    S_landing_reverse = S_A + S_F + S_FR + S_BR_reverse
    S_ground_roll_reverse = S_FR + S_BR_reverse
    
    # approach angle
    

    #Output in si units
    S_landing_si  = S_landing*0.3048
    S_ground_roll_si = S_ground_roll*0.3048
    S_landing_reverse_si  = S_landing_reverse*0.3048
    S_ground_roll_reverse_si = S_ground_roll_reverse*0.3048
    V_REF = V_REF*0.3048
    
    return(S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO, theta_app)


#Power should be imported, remove as input
#Cruise altitude can be set as a constant and removed as input
"landing weight and descent weight should be predefined"
"aproach angle should be predefined"
#T_static does not have to be calculated using Dprop etc...
#T_static might not have to be calculated at all

def main_updown(MTOW, S, A, e, CD_0, etha_p, etha_p_landing, P_TO, P, T_static, LD_max, hcruise, CL_0, CL_alpha, alpha_TO, CL_max_TO, CL_landing_max, CL_landing_td, D_prop, D_spinner, V_c, V_H):
    CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si = Takeoff(MTOW, S, A, e, CD_0, etha_p, P_TO, CL_0, CL_alpha, alpha_TO, CL_max_TO, D_prop, D_spinner, V_c, V_H)
    t_cruiseh, ROC_SL_si, h_absolute_si = Climb(MTOW, S, A, e, CD_0, etha_p, P, LD_max, hcruise)
    R_glide = Descent(MTOW, S, A, e, CD_0, LD_max, hcruise)
    S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO = Landing(MTOW, S, A, e, CD_0, etha_p_landing, CL_landing_max, T_static, theta_app, CL_landing_td) 
    return(CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si, t_cruiseh, ROC_SL_si, h_absolute_si, R_glide, S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO)

#0.965
"""
###############################################################################
############################  Printing  #######################################
###############################################################################

#PRINTING
#print "S_G =", S_G[loc_list], "ft"
#print "S_TR =", S_TR, "ft"
#print "S_C =", S_C, "ft"
#print "Take-off length =", S_TO_si , "m"
#print "Climb angle =", theta_climb, "degrees"
#print "Maximum rate of climb reachable at sea level =", ROC_SL,"m/s"
#print "Time to reach an altitude of", hcruise,"m is", t, "s"
#print "Gliding distance from cruise =", R_glide, "m, or", R_glide*0.000539956803, "nm"
#print "Ground roll distance =", S_ground_roll_si,"m"
#print "Total landing distance =", S_landing_si,"m"

#It should be noted that the maximum climb rate is calculated at MTOW with full fuel tanks


###############################################################################
###############################  PLOTS  #######################################
###############################################################################

#Thrust - airspeed plot
#plt.plot(V, T_V)
#plt.axis([0, 200, 0, 1500])
#plt.show()

#Speed -  take-off distance - lift-off speed plot
#plt.plot(S_G, V)
#plt.plot(S_G, VLOF_list)
#plt.axis([0, 1500, 0, 140])
#plt.show() 

#Altitude - ROC plot
#plt.plot(ROC_list, htab[0:len(ROC_list)])
#plt.show()

#ROC - Altitude plot
#plt.plot(htab[0:len(ROC_list)], ROC_list )
#plt.plot(htab[0:len(ROC_list)], polyfit)
#plt.show()
"""



