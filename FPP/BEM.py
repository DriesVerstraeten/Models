# -*- coding: utf-8 -*-
# FILE: BEM.py
"""
Created on Thu Jun 08 13:46:15 2017

@author: Simon

This file contains the classes used for the blade element theory and propulsion
calculations.

General assumptions
-------------------

The incoming velocity is always perpendicular to the propellor disk

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
        (in meters)
    
    chord: ndarray
        Contians the chord of the sections (in meters)
        
    beta: ndarray
        Contains the geometric pitch (beta) of the blade sections (in radians)
        
    thickness: ndarray
        Contains the thickness ratio of the sections (not scaled, no units)
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
    Class that defines airfoil properties
    
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
    
    def for_thickness(self, thickness):
        """
        Return interpolated lift & drag data for the given thickness.

        Parameters
        ----------
        thickness : float
            Fractional thickness

        """
        lift_drag = self.lift_drag_by_thickness(thickness)
        return lift_drag
        
def _strip_boundaries(radii):
    """
    Find the midpoints in the radii specified, the end points are appended to 
    the start and end.
    """
    radii = 1.0 * np.asarray(radii)
    midpoints = (radii[1:] + radii[:-1]) / 2
    return np.r_[radii[0], midpoints, radii[-1]]

def LTS(velocity, rps, radius):
    """
    Local total speed at a certain radius; 
    ||incoming vector + rotational vector||
    """
    return np.sqrt(velocity**2+(rps*radius*2*np.pi)**2)

def fast_interpolation(x, y, new_x):
    """
    from https://stackoverflow.com/a/13504757
    
    Computes an interpolation in multiple dimensions, but with greater 
    efficiency than multiple interp1d calls
    """
    from scipy.interpolate._fitpack import _bspleval
    f = interp1d(x, y, axis=-1, kind=3)
    xj,cvals,k = f._spline
    result = np.empty_like(new_x)
    for (i, j), value in np.ndenumerate(new_x):
        result[i, j] = _bspleval(value, x, cvals[:, i, j], k, 0)
    return result

def advance_angle(velocity, rps, radius):
    """
    Calculate the advance angle experienced by the rotorblades
    
    Input
    -----
    velocity:   float
        Free-stream velocity
    rps:        float
        Frequency of the rotor blades
    radius:     ndarray
        Radius at which the angle(s) is/are to be calculated
    """
    return np.arctan(velocity/(2*np.pi*rps*radius))
    

class BET(object):
    """
    Class that contains the Blade Element Theory model (WiP)
    
    Parameters
    ----------
    
    blade: blade object
        Contains the parameters of one blade
    
    spinner_r: float
        Defines the distance from the center of the disc to the root of the 
        blade or in other words, the spinner radius
        
    num_blades: int
        Contains the amount of blades in the propellor
        
    airfoil: airfoil object
        Contains the airfoil(s) polar data (AoA, cl, cd, cm)
    """
    
    def __init__(self, blade, spinner_r, num_blades, airfoil):
        self.blade = blade
        self.root_l = spinner_r
        self.num_blades = num_blades
        
        self.radii = spinner_r + np.asarray(self.blade.x)
        self.boundaries = _strip_boundaries(self.radii)
        
        #airfoil data
        self.alpha = airfoil.alpha
        self.lift_drag = np.array([
            airfoil.for_thickness(th)
            for th in self.blade.thickness])
        self._lift_drag_interp = fast_interpolation(
            airfoil.alpha, self.lift_drag, axis=1)
        self._last_factors = np.zeros((len(self.radii), 2)) #stores last coefficients
        
    def lift_drag(self, alpha):
        """
        Returns interpolated lift and drag values for the given angle of attack
        
        Input
        -----
        alpha:  array-like
            list of alpha values at each analysed section
        """
        alpha_vals = np.vstack((alpha, alpha)).T
        return self._lift_drag_interp(alpha_vals)
    
    def force_coeffs(self,advance_angle,pitch):
        """
        Returns the force coefficients out-of and in the plane (cT, cK)
        
        Parameters
        ----------
        inflow angle:   ndarray
            The velocity angle w.r.t. the disk
        pitch:          float
            The pitch setting of the prop
        """
        twist = self.blade.beta
        if len(twist) != len(advance_angle):
            raise ValueError("Shape mismatch, make sure the ")
        
        alpha = twist-advance_angle
        cl_cd = self.lift_drag(alpha)
        #construct array  for the conversion from lift to thrust
        cphi, sphi = np.cos(advance_angle), np.sin(advance_angle)
        A = np.array([[cphi, -sphi], [sphi, cphi]])
        return np.einsum("ijk,jk->ik", A, cl_cd)
        
    def forces(self,velocity, rho, rps, pitch):
        """
        Construct the force coefficients at each blade
        """
        r = self.radii
        chord = self.blade.chord
        lts = LTS(velocity, rps,r)
        phi = advance_angle(velocity, rps,r)
        force_coeffs = self.force_coeffs(phi,pitch)
        forces = 0.5*rho*lts**2*force_coeffs * chord
        return forces
        
    def thrust_torque(self,velocity, rps, rho, pitch = 0.0):
        """
        Generate the forces on the rotordisk
        """
        r = self.radii
        forces = self.forces(velocity, rho, rps)
        fx, fy = zip(*forces)
        thrust = self.num_blades * np.trapz(fx, x=r)
        torque = self.num_blades * np.trapz(-np.array(fy) * r, x=r)
        power = torque * rps * np.pi*2
        
        return thrust, torque, power
        
    
