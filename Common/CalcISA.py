# -*- coding: utf-8 -*-
"""
This script calculates the ISA conditions for a certain height,
currently editing it to take pressure or temperatures to calculate the height

Created on Tue Apr 28 09:09:24 2015

@author: Simon
"""

import numpy as np
#These are basic variables
P0 = 101325.0 #Pa
T0 = 288.15 #K
rho0 = 1.225 #kg/m^3
aTrop = -0.0065 #K/m
aStrat1 = 0.001
aStrat2 = 0.0028
aMeso1 = -0.0028
aMeso2 = -0.002
g0 = 9.80665 #m/s^2
R = 287.058 #J/(kg*K)
T0_v = 291.15 #K
C = 110.56 #K
mu_0 = 18.27*10**(-6) #Pa*s

#define the temperatures at the kinks in the temperature


#@11000 m
T11 = T0 + aTrop * 11000
PR11 = (T11/T0)**(-g0/(aTrop * R))
P11 = PR11*P0
rho11 = rho0 * PR11 * T0/T11

#@20000 m
T20 = T11
PR20 = np.exp(-g0/(R*T20)*9000)
P20 = P11*PR20
rho20 = rho11*PR20

#@32000 m
T32 = T20 + aStrat1 * 12000
PR32 = (T32/T20)**(-g0/(aStrat1 * R))
P32 = PR32*P20
rho32 = PR32*rho20*(T20/T32)

#@47000 m
T47 = T32 + aStrat2 * 15000
PR47 = (T47/T32)**(-g0/(aStrat2 * R))
P47 = PR47*P32
rho47 = PR47*rho32*(T32/T47)

#@51000 m
T51 = T47
PR51 = np.exp(-g0/(R*T51)*4000)
P51 = P47*PR51
rho51 = rho47*PR51

#@71000 m
T71 = T51 + aMeso1 * 20000
PR71 = (T71/T51)**(-g0/(aMeso1 * R))
P71 = P51*PR71
rho71 = PR71*rho51*(T51/T71)

#@84852 m
#T84 = T71 + aMeso2 * 13852

def ISA(h):
    #Define the temperature along the height
    #Troposphere
    if h<=11000 and h>=0:
        T1 = T0 + aTrop * h
        PR1 = (T1/T0)**(-g0/(aTrop * R))
        P1 = P0* PR1
        rho1 = rho0*PR1*T0/T1
        return np.array([T1, P1, rho1])
    
    #Tropopause
    elif h>11000 and h<= 20000:
        T1 = T11
        PR1 = np.exp(-g0/(R*T1)*(h-11000))
        P1 = P11 *PR1
        rho1 = rho11 *PR1
        return np.array([T1, P1, rho1])
    
    #Part one of the Stratosphere
    elif h>20000 and h<=32000:
        T1 = T20 + aStrat1 * (h-20000)
        PR1 = (T1/T20)**(-g0/(aStrat1 * R))
        P1 = P20* PR1
        rho1 = rho20*PR1*T20/T1
        return np.array([T1, P1, rho1])
    
    #Part two of the Stratosphere
    elif h>32000 and h<=47000:
        T1 = T32 + aStrat2 * (h-32000)
        PR1 = (T1/T32)**(-g0/(aStrat2 * R))
        P1 = P32* PR1
        rho1 = rho32*PR1*T32/T1
        return np.array([T1, P1, rho1])
    
    #Stratopause
    elif h>47000 and h<=51000:
        T1 = T47
        PR1 = np.exp(-g0/(R*T1)*(h-47000))
        P1 = P47 *PR1
        rho1 = rho47 *PR1
        return np.array([T1, P1, rho1])
        
    #Part one of the Mesosphere
    elif h>51000 and h<=71000:
        T1 = T51 + aMeso1 * (h-51000)
        PR1 = (T1/T51)**(-g0/(aMeso1 * R))
        P1 = P51* PR1
        rho1 = rho51*PR1*T51/T1
        return np.array([T1, P1, rho1])
    
    #Part two of the Mesosphere
    elif h>71000 and h<=84852:
        T1 = T71 + aMeso2 * (h-71000)
        PR1 = (T1/T71)**(-g0/(aMeso2 * R))
        P1 = P71* PR1
        rho1 = rho71*PR1*T71/T1
        return np.array([T1, P1, rho1])
        
    else:
        raise ValueError('{} does not lie between 0 and 84852 m, I cannot calculate ISA values there'.format(h))

def ISA_temp(h):
    '''
    Calculate the temperature at a certain height
    '''
    if h<=11000 and h>=0:
        T1 = T0 + aTrop * h
        return T1
    
    #Tropopause
    elif h>11000 and h<= 20000:
        return T11
    
    #Part one of the Stratosphere
    elif h>20000 and h<=32000:
        T1 = T0 + aStrat1 * h
        return T1
    
    #Part two of the Stratosphere
    elif h>32000 and h<=47000:
        T1 = T0 + aStrat2 * h
        return T1
    
    #Stratopause
    elif h>47000 and h<=51000:
        return T47
        
    #Part one of the Mesosphere
    elif h>51000 and h<=71000:
        T1 = T0 + aMeso1 * h
        return T1
    
    #Part two of the Mesosphere
    elif h>71000 and h<=84852:
        T1 = T0 + aMeso2 * h
        return T1

def ISA_press(h):
    T = ISA_temp(h)
    
    #Troposphere
    if h<=11000 and h>=0:
        P1 = P0*(T/T0)**(-g0/(aTrop * R))
        return P1
    
    #Tropopause
    elif h>11000 and h<= 20000:
        P1 = P11*np.exp(-g0/(R*T)*(h-11000))
        return P1
    
    #Part one of the Stratosphere
    elif h>20000 and h<=32000:
        P1 = P20*(T/T20)**(-g0/(aStrat1 * R))
        return P1
    
    #Part two of the Stratosphere
    elif h>32000 and h<=47000:
        P1 = P32*(T/T32)**(-g0/(aStrat2 * R))
        return P1
    
    #Stratopause
    elif h>47000 and h<=51000:
        P1 = P47*np.exp(-g0/(R*T)*(h-47000))
        return P1
        
    #Part one of the Mesosphere
    elif h>51000 and h<=71000:
        P1 = P51*(T/T51)**(-g0/(aMeso1 * R))
        return P1
    
    #Part two of the Mesosphere
    elif h>71000 and h<=84852:
        P1 = P71*(T/T71)**(-g0/(aMeso2 * R))
        return P1
    
def ISA_dens(h):
    T = ISA_temp(h)
    
    #Troposphere
    if h<=11000 and h>=0:
        rho1 = rho0*(T/T0)**(-g0/(aTrop * R)-1)
        return rho1
    
    #Tropopause
    elif h>11000 and h<= 20000:
        rho1 = P11*np.exp(-g0/(R*T)*(h-11000))
        return rho1
    
    #Part one of the Stratosphere
    elif h>20000 and h<=32000:
        rho1 = rho20*(T/T20)**(-g0/(aStrat1 * R)-1)
        return rho1
    
    #Part two of the Stratosphere
    elif h>32000 and h<=47000:
        rho1 = rho32*(T/T32)**(-g0/(aStrat2 * R)-1)
        return rho1
    
    #Stratopause
    elif h>47000 and h<=51000:
        rho1 = rho47*np.exp(-g0/(R*T)*(h-47000))
        return rho1
        
    #Part one of the Mesosphere
    elif h>51000 and h<=71000:
        rho1 = rho51*(T/T51)**(-g0/(aMeso1 * R)-1)
        return rho1
    
    #Part two of the Mesosphere
    elif h>71000 and h<=84852:
        rho1 = rho71*(T/T71)**(-g0/(aMeso2 * R)-1)
    
def dyn_visc(h):
    '''
    mu = mu_0*(T0+C)/(T+C)(T/T0)**(3/2)
    C = 120
    mu_o = 18.27*10e-6
    T0 = 291.15
    '''    
    T = ISA_temp(h)
    return mu_0*(T0_v+C)/(T+C)*((T/T0_v)**(3/2))

def speed_sound(h):
    T = ISA_temp(h)
    return np.sqrt(1.4*T*R)
    
