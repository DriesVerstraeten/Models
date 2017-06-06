#Importing modules
import numpy as np
from math import *

#Importing functions


#Field performance

#def TO_perf(mu_g, W_TO, S, D_p, P_TO, rho, CL_max_TO, V_c, V_H, etha_p, T_static, g):

#Initial test input values
CL_0 = 0.39
CL_alpha = 0.1
alpha_TO = 2 #deg
W_TO = 19620 #Newton (=2000 kg)
rho = 1.225 #kg/m3
S = 15 #m2
CL_max_TO = 1.1 
etha_p = 0.82
P_TO = 400 #BHP
V_c = 87.5 #m/s
V_H = 92.6 #m/s
T_static = 5827 #Newton
mu_g = 0.04 #Wet grass brakes off
g = 9.81 #m/s2
CD_min = 0.0350 #Cirrus SR22T value
A = 6
e = 0.8


#Function code
CL_TO = CL_0 + CL_alpha * alpha_TO
CD_TO = CD_min + (CL_TO**2/(pi*A*e))
V_S1 = sqrt((2*W_TO)/(rho*S*CL_max_TO))
V_LOF = 1.1 * V_S1
T_c = (etha_p*550*P_TO)/V_c
T_H = (etha_p*550*P_TO)/V_H   

a = np.matrix([[0, 0, 0, 1],
               [(V_c**3), (V_c**2), V_c, 1],
               [(3*V_c**2), (2*V_c), 1, 0],
               [(V_H**3), (V_H**2), V_H, 1]])
b = np.matrix([[T_static],
               [T_c],
               [-etha_p*325.8*(P_TO/(V_c**2))],
               [T_H]])
x = np.linalg.solve(a, b)
O = float(x[0])
B = float(x[1])
C = float(x[2])
E = float(x[3])



i = 0
time = 0
dt = 0.5
t = []
V = [0.]
L = []
D = []
q = []
a = [0.]
S1 = []
S2 = []
S_TO = [0.]
T_V = []

while i < 5:
    
    
    dynam_pressure = 0.5*rho*V[i]**2.   
    lift = dynam_pressure*S*CL_TO
    drag = dynam_pressure*S*CD_TO
    
    
    thrust = (O*V[i]**3.) + (B*V[i]**2.) + (C*V[i]) + E    
    acceleration = (thrust-drag-(mu_g*(W_TO - lift)))/(W_TO/g)
    airspeed = V[i-1] + a[i]*dt
    
    print V[i-1]
    print a[i]
    print airspeed
    print
    
    distance_1 = V[i-1] * dt
    distance_2 = 0.5 * a[i-1] * (dt**2)
    distance_total = S_TO[i-1] + distance_1 + distance_2
    
    
    q.append(dynam_pressure)
    L.append(lift)
    D.append(drag)
    a.append(acceleration)
    V.append(airspeed)
    S1.append(distance_1)
    S2.append(distance_2)
    S_TO.append(distance_total)
    t.append(time)
    T_V.append(thrust)        
    
    time = time + dt   
    i = i+1
        
        
#return(CL_TO, CD_TO, V_LOF)

#Climb performance


#Descent performance


