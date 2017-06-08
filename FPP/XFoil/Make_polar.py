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
    
    Inputs:
        RE:     float, Reynolds number
        M:      float, Mach number
    
    Outputs:
        RE_f:   str, Reynolds number
        M_f:    str, Mach number
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
    
    """    
    #Begin the command sequence
    issueCmd('load '+filename)  #load airfoil coordinates in XFoil
    issueCmd('GDES')            #enter GDES menu
    issueCmd('CADD')            #add points at corners
    issueCmd('')                #accept default input
    issueCmd('')                #accept default input
    issueCmd('')                #accept default input
    issueCmd('')                #return to top menu
    issueCmd('PPAR')            #enter paneling parameters menu
    issueCmd('N 299')           #change panel count to 299 (max number)
    issueCmd('')                #accept
    issueCmd('')                #accept
    issueCmd('')                #return to top level
    issueCmd('OPER')            #enter OPER menu
    issueCmd('VISC '+RE_s)      #change to viscous mode
    issueCmd('MACH '+M_s)       #change the mach number
    issueCmd('ITER 100')        #change iteration limit for polar accumulation
    issueCmd('PACC')            #change to polar mode
    issueCmd('../Polars/'+foilname +'_RE' + RE_f + '_M' + M_f +'.txt')   #specify filename
    issueCmd('')                #no dump needed
    issueCmd('ASEQ -10 10 0.5') #generate the linear, lower part
    issueCmd('ASEQ 10.1 25 0.1')#generate the upper part with the non-linear part
    issueCmd('')                #return to top level
    issueCmd('QUIT')            #quit XFOIL"""
    output = xfoil.communicate(string.join(['load '+filename, 'GDES', '', '', '', '', 
                                'PPAR', 'N200', '', '', '', 'OPER', 'VISC '+RE_s,
                                'MACH '+M_s, 'ITER 100', 'PACC', 
                                '../Polars/'+foilname +'_RE' + RE_f + '_M' + M_f +'.txt'
                                '', 'ASEQ -10 10 0.5', 'ASEQ 10.1 25 0.1',
                                'PACC'], '\n'))
    
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
    """
    RE_f, M_f = convert_float_to_str(RE,M)
    polar = load_foil_data(foilname, RE, M)
    alpha = polar[:,0]
    cl = polar[:,1]
    cd = polar[:,2]
    cm = polar[:,4]
    np.savez(r'../Polars/' + foilname +'_RE' + RE_f + '_M' + M_f +'.txt', AoA = alpha, cl = cl, cd = cd, cm = cm)