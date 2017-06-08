# -*- coding: utf-8 -*-
# FILE: BEM.py
"""
Created on Thu Jun 08 13:46:15 2017

@author: Simon
"""

import numpy as np

class blade():
    """
    Class that defines a blade
    
    Attributes
    ----------
    
    x: ndarray
        Contains the positions of the stations measured from blade root
    
    chord: ndarray
        Contians the chord of the sections
        
    beta: ndarray
        Contains the pitch of the blade sections
        
    thickness: ndarray
        Contains the thickness of the sections
    """
    def __init__(self, x, chord, beta, thickness):
        #making sure the attributes are ndarrays
        self.x = np.array(x)
        self.beta = np.array(beta)
        self.chord = np.array(chord)
        self.thickness = np.array(thickness)
        #check is all parameters are complete
        if not (len(x) == len(chord) == len(beta) == len(thickness)):
            raise ValueError("Shape mismatch")

class BET(object):
    """
    Class that contains the Blade Element Theory
    
    Attributes
    ----------
    
    
    """