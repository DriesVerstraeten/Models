#NOTES
#
#Power should still be implemented as a function of altitude
#Power used is mainly in BHP
#P_TO = P_BHP at sea level
#
#

#Importing modules
import numpy as np
from math import *
import matplotlib.pyplot as plt

#Importing functions
import Common.CalcISA as ISA

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

###############################################################################
##############################  INPUTS  #######################################
###############################################################################

#Initial test input values (IMPERIAL UNITS)

CL_0 = 0.39
CL_alpha = 0.1
alpha_TO = 2. #deg
W_TO = 3400. #lbf
rho = 0.002378 #slugs/ft3
S = 144.9 #ft2
CL_max_TO = 2.2 
etha_p = 0.85 #propeller efficiency
P_TO = 450. #BHP
V_c = 287. #ft/s - Cruising velocity
V_H = 314. #ft/s - Maximum velocity
D_prop = 5.577427822 #ft
D_spinner = 1.312335958 #ft

CD_0 = 0.0250 
A = 7 #Aspect Ratio
e = 0.85 #Oswald efficiency factor
LD_max = 16. #Maximum lift-over-drag


###############################################################################
#########################  Take-off performance  ##############################
###############################################################################


def TO_perf(CL_0, CL_alpha, alpha_TO, CD_0, A, e, W_TO, S, CL_max_TO, D_prop, D_spinner, V_c, V_H, etha_p, P_TO):

    #NOTES:
    #All input values and formulas are performed using imperial units
    #The method used is described in chapter 17 method #3
    
    # Not an input for the function, however an input locally
    rho_sealevel = 0.002378 #slugs/ft3
    mu_g_nobrakes = 0.04 #Ground friction coefficient
    g = 32.1740 #ft/s2
    
    #Converting the input values from si units to imperial units
    #W_TO = W_TO*0.224808943871 #Newton to lbf
    #S = S*10.7639104 #m2 to ft2
    #D_prop = D_prop*3.2808399 #m to ft
    #D_spinner = D_spinner*3.2808399 #m to ft
    #V_c = V_c*3.2808399 #m/s to ft/s
    #V_H = V_H*3.2808399 #m/s to ft/s
    #P_TO = P_TO*1.3410220888 #kW to BHP
    
    #Basic functions taken from Chapter 17 of the General Aviation book
    CL_TO = CL_0 + CL_alpha * alpha_TO #Can we get this as an output ?
    CD_TO = CD_0 + ((CL_TO**2)/(pi*A*e))
    V_S1 = sqrt((2*W_TO)/(rho_sealevel*S*CL_max_TO)) #ft/s
    V_LOF = 1.1 * V_S1 #ft/s
    T_c = (etha_p*550*P_TO)/V_c #lbf
    T_H = (etha_p*550*P_TO)/V_H #lbf  
    
    A_prop = 0.25*pi*D_prop**2
    A_spinner = 0.25*pi*D_spinner**2
    T_static = (0.85*((550*P_TO)**(2./3.))*((2.*rho_sealevel*A_prop)**(1./3.))*(1. - (A_spinner/A_prop))) #What power to use
    #T_static apart laten berekenen voor landing?
    #T_static als output van propulsion functie?
    
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
    dt = 0.1
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
    
    for i in range(0, 10000):
        if i == 0:
            time = 0
            thrust = E
            dynam_pressure = 0.5*rho_sealevel*V[i]**2
            lift = dynam_pressure*S*CL_TO
            drag = dynam_pressure*S*CD_TO
            acceleration = (thrust-drag-(mu_g_nobrakes*(W_TO - lift)))/(W_TO/g)
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
            thrust = (O*V[0]**3.) + (B*V[0]**2.) + (C*V[0]) + E
            dynam_pressure = 0.5*rho_sealevel*V[i-1]**2
            lift = dynam_pressure*S*CL_TO
            drag = dynam_pressure*S*CD_TO
            acceleration = (thrust-drag-(mu_g_nobrakes*(W_TO - lift)))/(W_TO/g)
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
            thrust = (O*V[i-1]**3.) + (B*V[i-1]**2.) + (C*V[i-1]) + E
            dynam_pressure = 0.5*rho_sealevel*V[i-1]**2
            lift = dynam_pressure*S*CL_TO
            drag = dynam_pressure*S*CD_TO
            acceleration = (thrust-drag-(mu_g_nobrakes*(W_TO - lift)))/(W_TO/g)
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
            
            
    VLOF_list = [V_LOF]*len(S_G) #Only useful for plotting
      
    loc = np.where(np.array(V) < V_LOF)
    loc_list = loc[0][-1]
    
    S_TR = (0.2156*((V_S1)**2)*((T_V[loc_list]/W_TO)-(1/(L[loc_list]/D[loc_list])))) #ft
    
    theta_climb = np.arcsin(((T_V[loc_list]/W_TO)-(1/(L[loc_list]/D[loc_list]))))*(180/pi) #degrees
    R_transition = (0.2156*(V_S1)**2) #ft
    h_transition = R_transition*(1. - np.cos(theta_climb*(pi/180))) #ft
    
    h_obst = 50 #ft, this is a CS23 requirement
    S_C = (h_obst-h_transition)/(np.tan(theta_climb*(pi/180))) #ft
    
    #SI units outputs
    S_TO_si = (S_G[loc_list] + S_TR + S_C)*0.3048
    V_LOF_si = V_LOF*0.3048
    V_S1_si = V_S1*0.3048

      
    return(CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si)


###############################################################################
########################  Climb performance  ##################################
###############################################################################

#def Climb(hcruise):

#Airspeed for best ROC
V_Y_list = []
hmax = 30000 #meter
h0 = 0
hcruise = 6000 #meter

for i in range(0, hmax):
    rho_h = Rhotab[i] * 0.00194032033 #Converted to slugs/ft3       
    V_Y = sqrt((2/rho_h)*(W_TO/S)*sqrt((1./(pi*A*e))/(3*CD_0)))
    V_Y_list.append(V_Y)

#Rate of climb
ROC_list = []

#If power changes according to height try changing the V_Y to a constant value

for i in range(0, len(V_Y_list)):
    ROC = 60*(((etha_p*550*P_TO)/W_TO)-(V_Y_list[i]*(1.1547/LD_max))) #ft/min, CHANGE P_TO to power dependent on altitude
    ROC_si = ROC*0.00508
    ROC_list.append(ROC_si)   
    ROC_SL = ROC_list[0]
    ROC_SL_si = ROC_SL*0.00508
    #if ROC<0.1 and ROC>-0.1:
        #print "Absolute ceiling =", htab[i], "m"
    #if ROC<100.1 and ROC>99.9:
        #print "Service ceiling =", htab[i], "m"
    
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

#return(t_cruiseh, ROC_SL)


###############################################################################
######################  Descent performance  ##################################
###############################################################################


#Best glide speed
V_BG_list = []

for i in range(0, hcruise):
    rho_h = Rhotab[i] * 0.00194032033 #Converted to slugs/ft3       
    V_BG = sqrt((2/rho_h)*(sqrt((1./(pi*A*e))/CD_0))*(W_TO/S)) #What weight are we going to take?
    V_BG_list.append(V_BG)

V_BG_list = V_BG_list[::-1] #Reversing the list, we're descending here

#Minimum glide angle
theta_min_glide = np.arctan(1./LD_max)


#Glide distance - provided the best glide speed and minimum glide angle
R_glide = hcruise * LD_max


###############################################################################
######################  Landing performance  ##################################
###############################################################################


#Adapt for landing at different altitudes

CL_landing_max = 2.2 
mu_g_brakes = 0.2 #Ground friction coefficient for wet grass while braking
T_landing = 0.07*T_static #lbf
T_landing_reverse_thrust = -0.6*T_static #lbf
theta_app = 3 #degrees
P_BHP_idle = 50 #TBD !!!!
CL_landing_td = 1.0 #CL after touchdown in landing config
CD_landing_td = 0.04 #CD after touchdown in landing config
etha_p_landing = 0.45

V_SO = sqrt((W_TO/S)*(2./rho)*(1./CL_landing_max)) #Landing stall speed in ft/s
V_REF = 1.3*V_SO #Approach speed
V_FLR = 1.3*V_SO #Flare speed in ft/s
V_TD = 1.1 * V_SO #Touchdown speed in ft/s
V_BR = 1.1 * V_SO #Brake speed in ft/s

L_landing_td = 0.5*rho*((V_BR/sqrt(2))**2)*S*CL_landing_td
D_landing_td = 0.5*rho*((V_BR/sqrt(2))**2)*S*CD_landing_td

h_f = 0.1512*(V_SO**2.)*(1.-np.cos(theta_app*(pi/180.)))
S_A = (h_obst - h_f)/(np.tan(theta_app*(pi/180.)))
S_F = 0.1512*(V_SO**2.)*np.sin(theta_app*(pi/180))
S_FR = V_TD
S_BR = abs(((V_BR**2.)*W_TO)/(2.*g*((sqrt(2.)*etha_p_landing*550.*(P_BHP_idle/V_BR))-D_landing_td-(mu_g_brakes*(W_TO-L_landing_td)))))


S_landing  = S_A + S_F + S_FR + S_BR
S_ground_roll = S_FR + S_BR

#Output in si units
S_landing_si  = S_landing*0.3048
S_ground_roll_si = S_ground_roll*0.3048


###############################################################################
############################  Printing  #######################################
###############################################################################

#PRINTING
print "S_G =", S_G[loc_list], "ft"
print "S_TR =", S_TR, "ft"
print "S_C =", S_C, "ft"
print "Take-off length =", S_TO_si , "m"
print "Climb angle =", theta_climb, "degrees"
print "Maximum rate of climb reachable at sea level =", ROC_SL,"m/s"
print "Time to reach an altitude of", hcruise,"m is", t, "s"
print "Gliding distance from cruise =", R_glide, "m, or", R_glide*0.000539956803, "nm"
print "Ground roll distance =", S_ground_roll_si,"m"
print "Total landing distance =", S_landing_si,"m"

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



