# -*- coding: utf-8 -*-
# FILE: BEM.py
"""
Created on Thu Jun 08 13:46:15 2017

@author: Simon

This file contains the classes used for the blade element theory and propulsion
calculations.

Modified from:  'https://github.com/ricklupton/py-bem'

General assumptions
-------------------

The incoming velocity is always perpendicular to the propellor disk

"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate._fitpack import _bspleval
import Common.CalcISA as ISA

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
        
    airfoils: dict of airfoil data arrays (AoA, cl, cd, cm)
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
        self.alpha = np.array(alpha)*np.pi/180
        
        lift_drag = np.dstack([
            [np.interp(alpha, data['alpha'], data['CL']) for data in airfoils],
            [np.interp(alpha, data['alpha'], data['CD']) for data in airfoils]
        ])
        self.lift_drag_by_thickness = interp1d(
            self.airfoils['thicknesses'], lift_drag, axis=0, copy=False)
    
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
'''
def fast_interpolation(x, y, new_x, axis):
    """
    from https://stackoverflow.com/a/13504757
    
    Computes an interpolation in multiple dimensions, but with greater 
    efficiency than multiple interp1d calls
    """
    from scipy.interpolate._fitpack import _bspleval
    f = interp1d(x, y, axis=axis, kind=3)
    xj,cvals,k = f._spline
    result = np.empty_like(new_x)
    for (i, j), value in np.ndenumerate(new_x):
        result[i, j] = _bspleval(value, x, cvals[:, i, j], k, 0)
    return result
'''
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

def _wrap_angle(theta):
    """
    Wraps the angle to [-pi, pi]
    """
    return (theta + np.pi) % (2 * np.pi) - np.pi

def cl_corr(cl, M):
    """
    Applies the Prandtl-Glauert correction to the Cl
    """
    return cl/np.sqrt(1-M)

def cd_corr(cd, M):
    """
    Applies the Frankl-Voishel correction to the Cd, and the drag divergence if
    M exceeds 0.7 (rough estimation)
    
    Del_mdd = 0.018
    """
    return cd*(0.000162*M**5 - 0.00383*M**4 + 0.0332*M**3
               - 0.118*M**2 + 0.0204*M + 0.996)

class fast_interpolation:
    def __init__(self, x, y, axis=-1):
        assert len(x) == y.shape[axis]
        self.x = x
        self.y = y
        self.axis = axis
        self._f = interp1d(x, y, axis=axis, kind='slinear', copy=False)

    def __getstate__(self):
        return dict(x=self.x, y=self.y, axis=self.axis)

    def __setstate__(self, state):
        self.x = state['x']
        self.y = state['y']
        self.axis = state['axis']
        self._f = interp1d(self.x, self.y, axis=self.axis,
                           kind='slinear', copy=False)

    def __call__(self, new_x):
        #assert new_x.shape == y.shape
        xj, cvals, k = self._f._spline.tck
        result = np.empty_like(new_x)
        for i, value in enumerate(new_x.flat):
            result.flat[i] = _bspleval(value, self.x, cvals[:, i], k, 0)
        return result
    


    

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
        
    height: float
        Contains the height at which the propellor is situated
    """
    
    def __init__(self, blade, spinner_r, num_blades, height, airfoil):
        self.blade = blade
        self.root_l = spinner_r
        self.num_blades = num_blades
        
        self.radii = spinner_r + np.asarray(self.blade.x)
        #self.boundaries = _strip_boundaries(self.radii)
        self.height = height
        
        #airfoil data
        self.alpha = airfoil.alpha
        self.lift_drag_dat = np.array([
            airfoil.for_thickness(th)
            for th in self.blade.thickness])
        self._lift_drag_interp = fast_interpolation(
            airfoil.alpha, self.lift_drag_dat, axis=1)
        
        
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
    
    def force_coeffs(self, inflow_ang, M_corr=False, LTS=None, T=None):
        """
        Returns the force coefficients out-of and in the plane (cT, cK)
        
        Parameters
        ----------
        inflow angle:   ndarray
            The velocity angle w.r.t. the disk
        """
        twist = self.blade.beta
        if len(twist) != len(inflow_ang):
            raise ValueError("Shape mismatch, make sure the amount of blade \
                             stations corresponds to the AoA")
        if M_corr:
            cl_cd = self.mach_factors(inflow_ang, LTS, T)
        else:
            cl_cd = self.lift_drag(inflow_ang)
        
        #construct array  for the conversion factors from lift to thrust
        cphi, sphi = np.cos(inflow_ang), np.sin(inflow_ang)
        A = np.array([[cphi, -sphi], [sphi, cphi]])
        return np.einsum("ijk,kj->ik", A, cl_cd)
        
    def forces(self, velocity, rho, rps, pitch = 0.0):
        """
        Construct the force coefficients at each blade
        """
        twist = self.blade.beta
        r = self.radii
        chord = self.blade.chord
        phi = advance_angle(velocity, rps,r)
        h = self.height
        
        alpha = twist-phi-pitch
        print 'alpha before: ', alpha
        ind_vel = np.ones(twist.shape)
        
        VE, alpha = self.induced_velocity(ind_vel, velocity, rps, alpha)
        
        print 'alpha after: ', alpha
        
        force_coeffs = self.force_coeffs(alpha, M_corr=True, LTS=VE, T=ISA.Temp(h))
        forces = 0.5*rho*VE**2*force_coeffs * chord
        return forces
        
    def thrust_torque(self, velocity, rps, rho, P, pitch = 0.0):
        """
        Generate the forces on the rotordisk
        """
        r = self.radii
        forces = self.forces(velocity, rho, rps, pitch)
        phi = advance_angle(velocity, rps,r)
        F = self.prandtl_corr(phi)
        thrust = self.num_blades * np.trapz(F * forces[0], x=r)
        torque = self.num_blades * np.trapz(F * -np.array(forces[1]) * r, x=r)
        power = np.abs(torque * rps * np.pi*2) 
        
        return thrust, torque, power
    
    def efficiency(self, velocity, rps, rho, pitch=0.0):
        """
        Calculates the efficiency of a propellor usinf the power coefficients
        """
        T, Q, P = self.thrust_torque(velocity, rps,rho, pitch)
        CT = T/((2*np.pi*rps)**2*rho*(2*self.radii[-1])**4)
        CP = P/((2*np.pi*rps)**3*rho*(2*self.radii[-1])**5)
        CQ = Q/((2*np.pi*rps)**2*rho*(2*self.radii[-1])**5)
        
        J = velocity/(2*np.pi*rps*2*self.radii[-1])
        
        eta = J*CT/CP
        return eta
    
    def induced_velocity(self, induced_vel, velocity, rps, alpha, max_iter=10):
        """
        Calculate the induced velocity component from an initial guess
        """
        running=True
        tol = 0.001
        alphai = 0.
        r = self.radii
        chord = self.blade.chord
        Nb = self.num_blades
        max_iter = 0
        T = ISA.Temp(self.height)
        while running:
            VE = LTS(velocity+induced_vel, rps, r)
            #print VE
            CL_CD = self.mach_factors(alpha-alphai, VE, T)
            
            
            #NEEDS CLEANING
            f = 8*np.pi*r/(chord*Nb)*induced_vel\
            - VE/(velocity + induced_vel)\
            * (CL_CD[:,0]*2*np.pi*rps*r - CL_CD[:,1]*(velocity + induced_vel))
            
            print 'grad: ', np.gradient(f)
            f_der = 8*np.pi*r/(chord*Nb)\
            - CL_CD[:,0]*2*np.pi*rps*r*(1/VE - VE/(velocity + induced_vel)**2)\
            * CL_CD[:,1]*(velocity + induced_vel)/VE
            
            wnew = induced_vel-f/f_der
            #print wnew
            #print induced_vel
            #print f/f_der
            diff = np.abs(wnew-induced_vel)
            #print diff
            max_iter += 1
            if np.all(diff)<=tol or max_iter==10:
                return LTS(velocity+wnew, rps, r), _wrap_angle(alphai)
            else:
                induced_vel = wnew
                alphai = np.arctan(wnew/VE)
                
    def mach_factors(self, alpha, LTS, T):
        """
        Defines the mach number and applies mach correction factors to the cl
        and cd.
        """
        cl_cd = self.lift_drag(alpha)
        Mach = LTS/ISA.speed_sound(T)
        #print Mach
        cl_mach = cl_corr(cl_cd[:,0],Mach)
        cd_mach = cd_corr(cl_cd[:,1],Mach)
        
        #print 'cl: ', cl_mach
        #print 'cd: ', cd_mach
        
        cl_cd[:,0] = np.where(Mach<0.3,cl_cd[:,0],cl_mach)
        cl_cd[:,1] = np.where(Mach<0.3,cl_cd[:,1],cd_mach)
        return cl_cd
    
    def prandtl_corr(self, phi):
        """
        Construcrs the prandtl factor 
        """
        Nb = self.num_blades
        r = self.radii
        if len(r) != len(phi):
            raise ValueError("Shape mismatch, make sure the amount of blade \
                             stations corresponds to the AoA")

        P_tip = Nb/2.*(r[-1]-r)/(r*np.sin(phi))
        P_root = Nb/2.*(r-r[0])/(r*np.sin(phi))
        F_tip = 2/np.pi*np.cos(np.exp(-1*P_tip))**(-1)
        F_root = 2/np.pi*np.cos(np.exp(-1*P_root))**(-1)
        return F_tip*F_root
        
