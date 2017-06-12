# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 15:07:52 2017

@author: Simon

Makes an XFOIL analysis and stores the polars in the polar folder

Naming convention: foilname_REXXXX_MXX.txt
    REXXXX is the reynolds number in tens of millions (e.g. RE0450 means a reynolds of 4.5 mil)
    MXX is the mach number (e.g. M03 means mach 0.3)
Polars should have a range of -10 to +25 AoA
"""

import psutil as sp
import numpy as np
from subprocess import PIPE
import string
import matplotlib.pyplot as plt


def start_XFOIL():
    global xfoil
    xfoil = sp.Popen(['xfoil.exe'],
      stdin=PIPE,
      stdout=PIPE,
      stderr=PIPE)
    return

def issueCmd(cmd,echo=True):
    xfoil.stdin.write(cmd+'\n')
    if echo:
        print cmd
        
def convert_float_to_str(RE,M):
    '''
    Convert the Reynolds number RE and Mach number M to floats suitable
    for our filenames.
    Reynolds number is converted to a three digit string, XXX, with the first
    digit equal to the tens of millions.
    Mach number is converted to a dwo digit string, XX, with the first digit 
    being the number before the point, e.g. 0.4 becomes '04'
    
    Input
    -----
    RE:     float
        Reynolds number
    M:      float
        Mach number
    
    Output
    ------
    RE_f:   str
        Reynolds number in string format
    M_f:    str
        Mach number in string format
    '''
    #Convert the Reynolds number to a string and get the tens of millions
    RE_f = str(np.round(RE))[:-7] #getting the parts higher than the thousands
    
    #checking whether we actually have the tens of millions
    while len(RE_f)<3:
        RE_f = '0' + RE_f
    
    M_f = str(np.round(M, decimals=1))
    
    return RE_f, M_f.replace('.','')
        
def generate_polar(foilname, RE, M, appendDAT = True):
    """
    Generate a polar with reynolds and mach number as specified
    """
    #generate filename
    
    RE_f, M_f = convert_float_to_str(RE,M)
    
    if appendDAT:
        filename = foilname + '.dat'
    else:
        filename = foilname
        foilname = foilname[:-4]  
       
    RE_s = str(np.round(RE/100000.)*100000)
    M_s = str(np.round(M, decimals=1))
    print '../Polars/'+foilname +'_RE' + RE_f + '_M' + M_f +'.txt'
    output = xfoil.communicate(string.join(['load '+filename, 'GDES', '', '', '', '', 
                                'PPAR', 'N200', '', '', '', 'OPER', 'VISC '+RE_s,
                                'MACH '+M_s, 'ITER 100', 'PACC', 
                                '../Polars/'+foilname +'_RE' + RE_f + '_M' + M_f +'.txt',
                                '', 'ASEQ -5 25 0.1',
                                'PACC', "", 'QUIT'], '\n'))
    
    return output

def run_analysis(foilname, RE, M):
    start_XFOIL()
    generate_polar(foilname, RE, M)
    return

def load_foil_data(foilname, RE, M):
    '''
    Loads a stored polar in XFOIL format, polars are found in the underlying
    polar map. First we construct the filename, including prepending the folder
    name and then we load the polar.
    Polars are in the following format:
        foilname_REXXXX_MXX.txt
    
    Input
    -----
    fname:  str
        filename of the foil requested without extension
    RE:     str
        Reynolds number at which the polar was constructed
    M:      float
        Mach number at which the polar was constructed
        
    Output
    ------
    data:   ndarray
        array containing the polar data
    '''
    #convert the RE and M to their resp strings
    RE_f, M_f = convert_float_to_str(RE,M)
    
    fname = r'../Polars/' + foilname +'_RE' + RE_f + '_M' + M_f +'.txt'
    data = np.genfromtxt(fname, dtype=float, skip_header=12)
    return data

def save_polar_npz(foilname, RE, M):
    """
    Save an XFoil polar to an NPZ file
    Alpha, CL, CD and CM are saved, their positions are:
    0, 1, 2, 4
    
    Input
    -----
    fname:  str
        filename of the foil requested without extension
    RE:     str
        Reynolds number at which the polar was constructed
    M:      float
        Mach number at which the polar was constructed
     
    """
    RE_f, M_f = convert_float_to_str(RE,M)
    polar = load_foil_data(foilname, RE, M)
    alpha = polar[:,0]
    cl = polar[:,1]
    cd = polar[:,2]
    cm = polar[:,4]
    np.savez(r'../Polars/' + foilname +'_RE' + RE_f + '_M' + M_f +'.txt', AoA = alpha, cl = cl, cd = cd, cm = cm)
    
def construct_thicknesses(foilname, thicknesses):
    """
    Constructs different thicknesses of a given airfoil by scaling the upper 
    and lower distribution over the camber line
    
    Input
    -----
    
    foilname:       str,
        filename of the aerofoil
    thicknesses:    ndarray
        Array of thicknesses to calculate the thickness distribution for
    """
    foil = load_foil(foilname)
    #get the upper and lower surfaces
    separation = np.where( np.logical_and(foil[:,0] == 0,foil[:,1] ==0))
    #print separation[0][0]
    s_upper = foil[:separation[0][0]+1][::-1]
    s_lower = foil[separation[0][0]:]
    
    #Get the camberline
    camber = (s_upper[:,1]-s_lower[:,1])/2+s_lower[:,1]
    #get the thickness distribution and scale it to the required relative thicknesses
    thick_distr = s_upper[:,1]-camber
    print 'Original thickness: ' ,np.amax(thick_distr)*2
    print 'New thicknesses: ', thicknesses
    new_foils_tc = np.zeros((len(thick_distr),len(thicknesses)))
    for i, thick in enumerate(thicknesses):
        new_foils_tc[:,i] = thick_distr / np.amax(thick_distr) * thick/2
    
    
    new_upper = new_foils_tc+np.transpose(np.array(camber, ndmin=2))
    new_lower = np.array(np.transpose(np.array(camber, ndmin=2))-new_foils_tc)[1:]
    
    new_foils = np.concatenate((new_upper[::-1,:],new_lower), axis=0)
    
    foil_out = np.concatenate((np.transpose(np.array(foil[:,0],ndmin=2)),new_foils,np.transpose(np.array(foil[:,1],ndmin=2))),axis=1)
    
    plt.plot(foil[:,0], foil[:,1], '-g')
    plt.plot(s_upper[:,0], camber, '-r')
    for i in range(len(foil_out[0,:])-2):
        plt.plot(foil_out[:,0], foil_out[:,i+1], '-b')
    
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()
    
    return foil_out, thicknesses

    
def load_foil(foilname, appendDAT=True):
    """
    Loads a coordinate file in selig format
    """
    if appendDAT:
        foilname = foilname + '.dat'
    foil = np.genfromtxt(foilname, skip_header=1)
    return foil