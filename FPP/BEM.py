# -*- coding: utf-8 -*-
# FILE: BEM.py
"""
Created on Thu Jun 08 13:46:15 2017

@author: Simon

This file contains the classes used for the blade element theory and propulsion
calculations.
"""

import numpy as np
from scipy.interpolate import interp1d

class blade(object):
    """
    Class that defines a blade
    
    Attributes
    ----------
    
    x: ndarray
        Contains the positions of the stations measured from blade root
    
    chord: ndarray
        Contians the chord of the sections
        
    beta: ndarray
        Contains the pitch (beta) of the blade sections in radians
        
    thickness: ndarray
        Contains the thickness of the sections
    """
    def __init__(self, x, chord, beta, thickness):
        #making sure the attributes are ndarrays
        self.x = np.array(x)
        self.beta = np.array(beta)
        self.chord = np.array(chord)
        self.thickness = np.array(thickness)
        #check if all parameters are complete and of the same shape
        if not (len(x) == len(chord) == len(beta) == len(thickness)):
            raise ValueError("Shape mismatch, make sure the inputs are of the same shape")
            
class airfoil(object):
    """
    Class that defines airfoil properties (WiP)
    
    It loads data from a '.npz' file with four arrays:
        
    airfoils: array of airfoil data arrays (AoA, cl, cd, cm)
    thickness: array of corresponding thickness values
    RE: Reynolds numbers for the corresponding airfoils
    M:  Mach number for the coresponding airfoils
    """
    def __init__(self, filename):
        self.filename = filename
        self.airfoils = np.load(filename)
        
        #make sure the alpha values are the same everywhere        
        alpha = []
        airfoils = self.airfoils['airfoils']
        for a in sorted(a for data in airfoils for a in data['alpha']):
            if alpha and abs(a - alpha[-1]) < 1e-5:
                continue
            alpha.append(a)
        self.alpha = np.array(alpha)
        
        lift_drag = np.dstack([
            [np.interp(alpha, data['alpha'], data['CL']) for data in airfoils],
            [np.interp(alpha, data['alpha'], data['CD']) for data in airfoils]
        ])
        self.lift_drag_by_thickness = interp1d(
            self.aerofoils['thicknesses'], lift_drag, axis=0, copy=False)
    
    

class BET(object):
    """
    Class that contains the Blade Element Theory model (WiP)
    
    Attributes
    ----------
    
    blade: blade object
        Contains the parameters of one blade
    
    spinner_r: float
        Defines the distance from the center of the disc to the root of the blade
        
    num_blades: int
        Contains the amount of blades in the propellor
        
    airfoil: airfoil object
        Contains the airfoil(s) polar data (AoA, cl, cd, cm)
    """
    
    def __init__(self, blade, spinner_r, num_blades, airfoil):
        self.blade = blade
        self.root_l = spinner_r
    
