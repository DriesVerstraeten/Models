# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 17:25:47 2017

@author: Simon

This is the propulsion file for the model

Inputs from chosen engine: AI-450S or M250 B17F (for now) (450 shp)

Other candidates: 
    GE H80 (or H75) (high performance)
    Some diesel (Continental CD300 or so) (economical option, no aerobatics)
    PT6A (medium version)
    
Using Blade Element Theory, I calculate the propellor chacteristics.
    (mainly efficiency and power available)
    
Prop geometry is saved in a blade class inside BEM.py
BEM.py contains all the relevant classes nd functions to perform the blade 
element (momentum?) theory
"""


import numpy as np
import Common.CalcISA as ISA
import BEM
import Init_Parameters as p
from scipy.interpolate import interp1d
import os

#could be stored in the parameter file
D = 1.8         #Prop diameter in meter
N = 4           #Number of blades
d_spin = 0.4    #Spinner diameter

#initial parameters
dx = 0.1        #Element length in meter


def chord_distr(x, D):
    """
    Defines a chord distribution over a vector x
    """
    unit_chord = np.array([.1, .14, .16, .16])
    x1 = np.array([0,.33,.66,1])
    
    chord_distr = interp1d(x1, unit_chord)
    x_scaled = x/x[-1]
    chord = chord_distr(x_scaled)*D
    return chord

def thick_distr(x):
    """
    Defines a thickness distribution over a vector x
    """
    unit_thick = np.array([0.2, 0.18, 0.16, 0.14, 0.12, 0.10, 0.10])
    x1 = np.array([0, 0.1, 0.2, 0.4, 0.6, 0.8, 1])
    thick_distr = interp1d(x1, unit_thick)
    x_scaled = x/x[-1]
    thick = thick_distr(x_scaled)
    return thick


def Analyse_prop(airfoil_path, h, V, rps, pitch = 0.0):
    """
    Analyse a propellor
    
    Input
    -----
    airfoil_path:   str
        Path to airfoil used in the analysis
    h:              float
        Height at which the analysis
    V:              float
        Velocity at which to calculate the analysis
    rps:            float
        Revolutions per second of the prop
    pitch:          float
        Additional pitch
    """
    x = np.arange(0,(D-d_spin)/2+dx,dx) #blade length
    chord_ref = interp1d(x,chord_distr(x, x[-1]))
    
    chord = chord_ref(x)
    beta = np.linspace(1.2,0.4,len(x))
    thickness = thick_distr(x)
    
    prop_blade = BEM.blade(x, chord, beta, thickness)
    
    e851 = BEM.airfoil(airfoil_path)
    
    BEM_analysis = BEM.BET(prop_blade, d_spin/2, N, e851)
    rho = ISA.Dens(h)
    thrust, torque, power = BEM_analysis.thrust_torque(V, rps, rho, pitch)
    eta_p = 2/(1+np.sqrt(1+thrust/(0.5*rho*V**2*D**2)))
    
    power_available = thrust*V
    eta_prop = power_available/power
    
    print eta_p
    print eta_prop
    return thrust, torque, power_available
