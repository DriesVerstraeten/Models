# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 15:20:45 2017

@author: Maxime
"""

def calcOEW(W_fuselage,W_wing,W_tail, mtow,Mfuelmax): # input initial value for mtow
    W_avionics = 50
    W_engine = 200
    
    OEW = W_fuselage + W_wing + W_tail  +  W_avionics + W_engine
    
    W_gear, W_Fuel_sys, W_flightcontrol, W_els, W_api, W_fur = Class2est_simplified(mtow,Mfuelmax,W_avionics)
    OEW += W_gear + W_Fuel_sys + W_flightcontrol + W_els + W_api + W_fur
    return OEW # [kg]



def Class2est_simplified(mtow,Mfuelmax,W_avionics):

    N_pax = 4 # People capacity
    M_D = 0.4 # Mach number

    sep_tanks = 2 # Number of fuel tanks
    Ne = 1 # Number of engines         
    
    Nult = 9 # Ultimate load factor
    lsm = 0.75 # Length of landing gear strut [m]
    integralfrac = 1 # Fuel tank: Qtot/(Qtot+Qint)

    W_gear =(0.054*((lsm/.3048)**0.501)*(((mtow*2.204)*Nult)**0.684))*0.45359237
    W_Fuel_sys =(2.49*((Mfuelmax*2.204/6.55)**0.6*(1./(1+integralfrac))**0.3*(sep_tanks)**0.2*(Ne)**0.13)**1.21)*0.45359237
    W_flightcontrol =1.066*((mtow*2.204)**0.626)*0.453592
    W_els = 426.0*((((W_Fuel_sys/0.45359237)+(W_avionics/0.45359237))/1000)**0.51)*0.45359237
    W_api = 0.265*((mtow*2.204)**0.52)*(N_pax**0.68)*((W_avionics/0.453592)**0.17)*(M_D**0.08)*0.453592
    W_fur = 0.412*((N_pax)**1.145)*(((mtow*2.204))**0.489)*0.453592
 
    return W_gear, W_Fuel_sys, W_flightcontrol, W_els, W_api, W_fur


mtow = 1700
Mfuelmax = 450

W_fuselage = 200
W_wing = 150
W_tail = 100

OEW = calcOEW(W_fuselage,W_wing,W_tail, mtow,Mfuelmax)

W_avionics = 450
W_gear, W_Fuel_sys, W_flightcontrol, W_els, W_api, W_fur = Class2est_simplified(mtow,Mfuelmax,W_avionics)

print OEW