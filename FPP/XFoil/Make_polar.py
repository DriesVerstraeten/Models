# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 15:07:52 2017

@author: Simon

Makes an XFOIL analysis and stores the polars in the polar folder

Naming convention: foilname_REXXXX_MXX.txt
    REXXXX is the reynolds number in tens of millions 
    (e.g. RE0450 means a reynolds nr of 4.5 mil)
    MXX is the mach number (e.g. M03 means mach 0.3)
Polars should have a range of -10 to +20 AoA
"""

import psutil as sp
import numpy as np
from subprocess import PIPE
import string
import matplotlib.pyplot as plt
import numpy.core.defchararray as defchar
from scipy.interpolate import interp1d

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
    print 'Saving to: ' + '../Polars/'+foilname +'_RE' + RE_f + '_M' + M_f +'.txt'
    output = xfoil.communicate(string.join(['load '+filename, 
                                            'GDES', '', '', '', '', 
                                            'PPAR', 'N 200', '', '',
                                            '', 'OPER', 'VISC '+RE_s,
                                            'MACH '+M_s, 
                                            'ITER 200',
                                            
                                            'PACC', 
                                '../Polars/'+foilname +'_RE' + RE_f + '_M' + M_f +'.txt',
                                '', 'ASEQ -5 20 0.1',
                                'PACC',
                                'INIT',
                                "", 'QUIT'], '\n'))
    
    return output

def run_analysis(foilname, RE, M):
    start_XFOIL()
    generate_polar(foilname, RE, M)
    return

def load_foil_data(foilname, RE_f, M_f):
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
    M:      str
        Mach number at which the polar was constructed
        
    Output
    ------
    data:   ndarray
        array containing the polar data
    '''
    fname = r'../Polars/' + foilname +'_RE' + RE_f + '_M' + M_f +'.txt'
    data = np.genfromtxt(fname, dtype=float, skip_header=12)
    return data

def save_1polar_npz(foilname, RE, M):
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
    polar = load_foil_data(foilname, RE_f, M_f)
    print len(polar)
    if len(polar)<1:
        np.savez(r'../Polars/' + foilname +'_RE' + RE_f + '_M' + M_f)
    else:
        alpha = polar[:,0]
        cl = polar[:,1]
        cd = polar[:,2]
        cm = polar[:,4]
        np.savez(r'../Polars/' + foilname +'_RE' + RE_f + '_M' + M_f, 
                 alpha = alpha, CL = cl, CD = cd, CM = cm)
    
def construct_thicknesses(foilname, thicknesses, plot=False):
    """
    Constructs different thicknesses of a given airfoil by scaling the upper 
    and lower distribution over the camber line
    
    Input
    -----
    foilname:       str,
        filename of the aerofoil
    thicknesses:    ndarray
        Array of thicknesses to calculate the thickness distribution for
        
    Output
    ------
    foil_out:       narray
        Two-dimensional array that contains the x-locations on the first
        collumn, the different thicknesses on the next collumns and the 
        original on the last collumn
    thicks_out:     ndarray
        Array containing the different thickensses at whick the airfoil is 
        reconstructed
    """
    foil = load_foil(foilname)
    #get the upper and lower surfaces
    separation = np.where( np.logical_and(foil[:,0] == 0,foil[:,1] ==0))
    #print separation[0][0]
    s_upper = foil[:separation[0][0]+1][::-1]
    s_lower = foil[separation[0][0]:]
    print foil.shape
    
    #Get the camberline
    if len(s_upper[:,0])==len(s_lower[:,0]):
        camber = (s_upper[:,1]-s_lower[:,1])/2+s_lower[:,1]
        #get the thickness distribution
        thick_distr = s_upper[:,1]-camber
    else:
        #some more work is needed to get the coordinates in the right spacing
        #first interpolate both upper and lower side to specified densities
        x1 = np.arange(0,0.01,0.0005)
        x2 = np.arange(0.01,0.05,0.01)
        x3 = np.arange(0.05,0.95,0.025)
        x4 = np.arange(0.95,1.001,0.01)
        x = np.concatenate((x1, x2, x3, x4),axis=0)
        s_upper_interp = interp1d(s_upper[:,0],s_upper[:,1])
        s_upper = np.stack((x,s_upper_interp(x)),axis=-1)
        s_lower_interp = interp1d(s_lower[:,0],s_lower[:,1])
        s_lower = np.stack((x,s_lower_interp(x)),axis=-1)
        
       
        
        #reconstruct the foil
        foil = np.concatenate((s_upper[::-1,:],s_lower[1:]), axis=0)
        #now we can calculate the camber
        camber = (s_upper[:,1]-s_lower[:,1])/2+s_lower[:,1]
        #get the thickness distribution
        thick_distr = s_upper[:,1]-camber
        print foil.shape
        print x.shape
    print 'Original thickness: ' ,np.amax(thick_distr)*2
    print 'New thicknesses: ', thicknesses
    #scale the thickness to the requirested thicknesses
    new_foils_tc = np.zeros((len(thick_distr),len(thicknesses)))
    for i, thick in enumerate(thicknesses):
        new_foils_tc[:,i] = thick_distr / np.amax(thick_distr) * thick/2
    
    
    new_upper = new_foils_tc+np.transpose(np.array(camber, ndmin=2))
    new_lower = np.array(np.transpose(np.array(camber, ndmin=2))-new_foils_tc)[1:]
    
    new_foils = np.concatenate((new_upper[::-1,:],new_lower), axis=0)
    
    foil_out = np.concatenate((np.transpose(np.array(foil[:,0],ndmin=2)),
                               new_foils,np.transpose(np.array(foil[:,1],
                                                               ndmin=2))),
                                                        axis=1)
    if plot:
        plt.plot(foil[:,0], foil[:,1], '-g')
        plt.plot(s_upper[:,0], camber, '-r')
        for i in range(len(foil_out[0,:])-2):
            plt.plot(foil_out[:,0], foil_out[:,i+1], '-b')
        
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()
    
    return foil_out, np.concatenate((thicknesses,np.array(np.amax(thick_distr)*2,
                                                          ndmin=1)),axis=0)

def load_foil(foilname, appendDAT=True):
    """
    Loads a coordinate file in selig format
    """
    if appendDAT:
        foilname = foilname + '.dat'
    foil = np.genfromtxt(foilname, skip_header=1)
    return foil

def save_polar_set(foilname,RE,M):
    """
    Runner file that saves different analyses of polars of different 
    thicknes variations of a single airfoil to a combined '.npz' file
    
    Input
    -----
    foilname:   str
        Filename of the airfoil to analyse
    RE:         ndarray
        One-dimensional array containing the list of Reynolds numbers to 
        analyse
    M:          ndarray
        One-dimensional array containing the list of Mach numbers to analyse
    """    
    #construct different thicknesses around the camber line
    thicks = np.arange(0.10,0.40, 0.02)
    new_foils, thick_out = construct_thicknesses(foilname, thicks, plot=False)
    thick_out = defchar.replace(thick_out.astype(str),'.','')
    thick_out = defchar.replace(thick_out,'0','', count=1)
    
    #Analyse the constructed foils over a Reynolds and mach range,
    #construct the filenames and resave them as .npz files
    fnames = []
    data_foil = []
    data_thickness = []
    data_RE = []
    data_M = []
    
    for i in range(len(new_foils[0,:])-1):
        save_foils = np.stack((new_foils[:,0].reshape(1,len(new_foils[:,0])),
                               new_foils[:,i+1].reshape(1,len(new_foils[:,0]))),
        axis=-1)
        np.savetxt(foilname+thick_out[i] + '.dat',save_foils[0],
                   header=foilname+thick_out[i],comments='', fmt = '%.4f')
        
        for j in RE:
            for k in M:
                RE_f, M_f = convert_float_to_str(j,k)
                fnames.append(foilname+thick_out[i] +'_RE' + RE_f + '_M' + M_f)
                run_analysis(foilname+thick_out[i],j,k)
                save_1polar_npz(foilname+thick_out[i],j,k)
                polar = np.load(r'../Polars/'+foilname+thick_out[i]+'_RE' +
                                RE_f + '_M' + M_f + '.npz')
                data_foil.append(polar)
                data_thickness.append(thick_out[i])
                data_RE.append(j)
                data_M.append(k)
                
    np.savez(foilname,airfoils = data_foil, thickness = data_thickness,
             RE = data_RE, M = data_M)
    
    return

def save_existing_polar(fnames, savename, thicks):
    """
    Save an existing polar set to a '.npz' datafile
    
    Input
    -----
    fnames:     array-like
        List or array of strings that contains all the filenames
    """
    fnames = np.asarray(fnames)
    print fnames.dtype
    print 
    fnames = defchar.add('../Polars/',fnames)
    savedata = []
    for fname in fnames:
        polar = np.genfromtxt(fname, skip_header=11)
        alpha = polar[:,0]
        cl = polar[:,1]
        cd = polar[:,2]
        cm = polar[:,4]
        savedata.append(dict(alpha = alpha, CL = cl, CD = cd, CM = cm))
    print savedata
    np.savez('../Polars/'+savename, airfoils = savedata, thicknesses = thicks)
        