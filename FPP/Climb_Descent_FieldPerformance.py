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

################################################################################


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
V_c = 287. #ft/s
V_H = 314. #ft/s
D_prop = 5.577427822 #ft
D_spinner = 1.312335958 #ft
mu_g = 0.04 #Ground friction coefficient
g = 32.1740 #ft/s2
CD_0 = 0.0250 
A = 7 #Aspect Ratio
e = 0.85 #Oswald efficiency factor


LD_max = 16. #Maximum lift-over-drag

################################################################################
#Field performance

#def TO_perf(mu_g, W_TO, S, D_p, P_TO, rho, CL_max_TO, V_c, V_H, etha_p, T_static, g):

#NOTES:
#All input values and formulas are performed using imperial units
#The method used is described in chapter 17 method #3
#

#Function code

#Basic functions taken from Chapter 17 of the General Aviation book
CL_TO = CL_0 + CL_alpha * alpha_TO
CD_TO = CD_0 + (CL_TO**2/(pi*A*e))
V_S1 = sqrt((2*W_TO)/(rho*S*CL_max_TO))
V_LOF = 1.1 * V_S1
T_c = (etha_p*550*P_TO)/V_c
T_H = (etha_p*550*P_TO)/V_H   

A_prop = 0.25*pi*D_prop**2
A_spinner = 0.25*pi*D_spinner**2
T_static = (0.85*((550*P_TO)**(2./3.))*((2.*rho*A_prop)**(1./3.))*(1. - (A_spinner/A_prop)))

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
        dynam_pressure = 0.5*rho*V[i]**2
        lift = dynam_pressure*S*CL_TO
        drag = dynam_pressure*S*CD_TO
        acceleration = (thrust-drag-(mu_g*(W_TO - lift)))/(W_TO/g)
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
        dynam_pressure = 0.5*rho*V[i-1]**2
        lift = dynam_pressure*S*CL_TO
        drag = dynam_pressure*S*CD_TO
        acceleration = (thrust-drag-(mu_g*(W_TO - lift)))/(W_TO/g)
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
        dynam_pressure = 0.5*rho*V[i-1]**2
        lift = dynam_pressure*S*CL_TO
        drag = dynam_pressure*S*CD_TO
        acceleration = (thrust-drag-(mu_g*(W_TO - lift)))/(W_TO/g)
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
        
        
plt.plot(V, T_V)
plt.axis([0, 200, 0, 1500])
plt.show()

VLOF_list = [V_LOF]*len(S_G)

plt.plot(S_G, V)
plt.plot(S_G, VLOF_list)
plt.axis([0, 1500, 0, 140])
plt.show()        

loc = np.where(np.array(V) < V_LOF)
loc_list = loc[0][-1]

S_TR = (0.2156*((V_S1)**2)*((T_V[loc_list]/W_TO)-(1/(L[loc_list]/D[loc_list])))) #ft

theta_climb = np.arcsin(((T_V[loc_list]/W_TO)-(1/(L[loc_list]/D[loc_list]))))*(180/pi) #degrees
R_transition = (0.2156*(V_S1)**2) #ft
h_transition = R_transition*(1. - np.cos(theta_climb*(pi/180))) #ft

h_obst = 40 #ft
S_C = (h_obst-h_transition)/(np.tan(theta_climb*(pi/180))) #ft



#SI units outputs
S_TO_si = (S_G[loc_list] + S_TR + S_C)*0.3048
V_LOF_si = V_LOF*0.3048

#PRINTING
print "S_G =", S_G[loc_list], "ft"
print "S_TR =", S_TR, "ft"
print "S_C =", S_C, "ft"
print "Take-off length =", (S_G[loc_list] + S_TR + S_C)*0.3048 , "m"
print      
print "Climb angle =", theta_climb, "degrees"
print "___________________________________________________________________"
        
#return(CL_TO, CD_TO, V_LOF)


################################################################################


#Climb performance


#Airspeed for best ROC
V_Y_list = []
hmax = 20000 #meter
h0 = 0
hcruise = 6000 #meter

for i in range(0, hmax):
    rho_h = Rhotab[i] * 0.00194032033 #Converted to slugs/ft3       
    V_Y = sqrt((2/rho_h)*(W_TO/S)*sqrt((1./(pi*A*e))/(3*CD_0)))
    V_Y_list.append(V_Y)

#Rate of climb
ROC_list = []

for i in range(0, len(V_Y_list)):
    ROC = 60*(((etha_p*550*P_TO)/W_TO)-(V_Y_list[i]*(1.1547/LD_max))) #ft/min, CHANGE P_TO to power dependent on altitude
    ROC_si = ROC*0.00508
    ROC_list.append(ROC_si)   
    ROC_SL = ROC_list[0]
    ROC_SL_si = ROC_SL*0.00508

print
print "Maximum rate of climb reachable at sea level =", ROC_SL,"m/s"

plt.plot(ROC_list, htab[0:len(ROC_list)])
plt.show()
#plt.axis([14,18,0,7000])

plt.plot(htab[0:len(ROC_list)], ROC_list )

bb = ROC_list[0]
aa = (ROC_list[-1]-ROC_list[0])/(htab[len(ROC_list)]-htab[0])

       
#Time to altitiude & service ceiling
ROC_a = (aa*(hcruise-h0))/(np.log(aa*hcruise + bb) - np.log(aa*h0 + bb))
t = (hmax - h0)/ROC_a
print
print "Time to cruise altitude =", t, "s"



#Descent performance














