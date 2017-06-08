#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:17:57 2017

@author: driesverstraeten
"""

import Init_Parameters as p
import Wing_model as wm
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt


t = p.t_skin


def wingbox_MOI(airfoil_file_name):
    airfoil_coordinates = np.genfromtxt(airfoil_file_name,skip_header=1)
    
    ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
    #INTERPOLATION OF AIRFOIL DATAPOINTS
    
    xcoordinates = np.zeros(len(airfoil_coordinates))
    ycoordinates = np.zeros(len(airfoil_coordinates))
    
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1]
    
    fit1 = np.polyfit(xcoordinates[0:134],ycoordinates[0:134],25)
    f1 = np.poly1d(fit1)
    
    fit2 = np.polyfit(xcoordinates[133:256],ycoordinates[133:256],25)
    f2 = np.poly1d(fit2)
    
    #arc length calculator: http://www.emathhelp.net/calculators/calculus-2/definite-integral-calculator/?a=3/20&b=23/40&f=sqrt%28%286958%2Ax%20-%203239%29%5E2%20%2B%20100000000%29/10000&steps=on&var=x
    
    front_spar = 0.15
    back_spar = 0.575
    
    arc_length_US = 0.427622 #upper skin
    arc_length_LS = 0.425 #lower skin
    
    front_spar_length = np.abs(f1(front_spar)) + np.abs(f2(front_spar))
    back_spar_length = np.abs(f1(back_spar)) + np.abs(f2(back_spar))
    
    front_spar_area = np.abs(f1(front_spar)) + np.abs(f2(front_spar))*t
    back_spar_area = np.abs(f1(back_spar)) + np.abs(f2(back_spar))*t
    
    total_area = (arc_length_US + arc_length_LS)*t + front_spar_area + back_spar_area
                
                 
    sections = 1000.             
    dx = 1./sections             
    area_section = dx * t
    x = np.arange(front_spar,back_spar+dx,dx)
    
    y_NA = (area_section * (sum(f1(x)) + sum(f2(x))) + front_spar_area*(f1(front_spar) + f2(front_spar))/2 + back_spar_area*(f1(back_spar) + f2(back_spar))/2)/ total_area 
    x_NA = (area_section * sum(x) * 2 + front_spar_area * front_spar + back_spar_area*back_spar) / total_area
                 
    Ixx = 2*dx*t**3/12 + area_section * sum(f1(x)**2 + f2(x)**2) + t/12 * (front_spar_length**3 + front_spar_area*(f1(front_spar) + f2(front_spar - y_NA))**2 + back_spar_length**3 + back_spar_area*(f1(back_spar) + f2(back_spar - y_NA))**2)
    Iyy = (front_spar_area * (front_spar-x_NA)**2) + (back_spar_area * (back_spar-x_NA)**2) + area_section*(sum(x**2)-2*x_NA*sum(x)+len(x)*x_NA**2)
    
    return Ixx, Iyy

wingbox_MOI('foil1_modified.dat')
      
#plt.plot(x,f1(x))
#plt.plot(x,f2(x))
#plt.show()


"""

########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#WINGBOX SPAR LENGTHS

def wingbox_spar():
    
    hfront = []
    hback = []
    hfront_per = 0.15
    hback_per = 0.57
    y_front_coordinates = []
    y_back_coordinates = []
    
    airfoil_coordinates = np.round(np.genfromtxt('foil1_modified.dat',skip_header=1),2)
    
    for i in range(len(airfoil_coordinates)):
        if airfoil_coordinates[i][0] == hfront_per:
            if len(hfront) == 2:
                continue
            hfront.append(np.abs(airfoil_coordinates[i][1]))
            h_front = sum(hfront)
            y_front_coordinates.append(airfoil_coordinates[i][1])
            
        if airfoil_coordinates[i][0] == hback_per:
            if len(hback) == 2:
                continue
            hback.append(np.abs(airfoil_coordinates[i][1]))
            h_back = sum(hback)
            y_back_coordinates.append(airfoil_coordinates[i][1])

        wingbox_width = hback_per - hfront_per

    return h_front, h_back, y_front_coordinates, y_back_coordinates, wingbox_width, hfront_per, hback_per

#print wingbox_spar()[0:4]


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#NEUTRAL AXIS LOCATION

def wingbox_NA():

    #DATUMS
    y_datum = -0.07
    x_datum = 0.
    
    #DIMENSIONS
    bottom_skin_length = np.sqrt((wingbox_spar()[3][1]-wingbox_spar()[2][1])**2 + wingbox_spar()[4]**2)
    top_skin_length = np.sqrt((wingbox_spar()[3][0]-wingbox_spar()[2][0])**2 + wingbox_spar()[4]**2)
    
    #print bottom_skin_length, top_skin_length
    
    #Y - CG LOCATIONS
    front_spar_ycg = wingbox_spar()[0]/2 + np.abs(y_datum - wingbox_spar()[2][1])
    back_spar_ycg = wingbox_spar()[1]/2 + np.abs(y_datum - wingbox_spar()[3][1])
    top_skin_ycg = (wingbox_spar()[2][0]+wingbox_spar()[3][0])/2+np.abs(y_datum)
    bottom_skin_ycg = (wingbox_spar()[2][1]+wingbox_spar()[3][1])/2+np.abs(y_datum)
    
    #print front_spar_ycg, back_spar_ycg, top_skin_ycg, bottom_skin_ycg
    
    
    #X - CG LOCATIONS
    front_spar_xcg = wingbox_spar()[5]
    back_spar_xcg = wingbox_spar()[6]
    top_skin_xcg = wingbox_spar()[6]-wingbox_spar()[5]
    bottom_skin_xcg = wingbox_spar()[6]-wingbox_spar()[5]
    
    #print front_spar_xcg, back_spar_xcg, top_skin_xcg, bottom_skin_xcg

    #AREAS
    front_spar_area = wingbox_spar()[0]*t
    back_spar_area = wingbox_spar()[1]*t
    top_skin_area = top_skin_length*t  
    bottom_skin_area = bottom_skin_length*t 
    
    #print front_spar_area, back_spar_area, top_skin_area, bottom_skin_area
    
    total_area = front_spar_area + back_spar_area + top_skin_area + bottom_skin_area                          
    
    #NEUTRAL AXIS
    y_NA = (front_spar_area*front_spar_ycg + back_spar_area*back_spar_ycg + top_skin_area*top_skin_ycg + bottom_skin_area*bottom_skin_ycg)/total_area + y_datum 
    x_NA = (front_spar_area*front_spar_xcg + back_spar_area*back_spar_xcg + top_skin_area*top_skin_xcg + bottom_skin_area*bottom_skin_xcg)/total_area       
           
    return y_NA, x_NA, front_spar_area, front_spar_ycg, back_spar_area, back_spar_ycg, top_skin_length, bottom_skin_length
    #return {'y_NA':y_NA, 'x_NA':x_NA, front_spar_area, front_spar_ycg, back_spar_area, back_spar_ycg, top_skin_length, bottom_skin_length}


def wingbox_MOI():
    
    sections = 10000.
    ds = 1./sections
    s1 = np.arange(0,wingbox_NA()[6]+ds,ds) #top
    s2 = np.arange(0,wingbox_NA()[7]+ds,ds) #bottom 
    alpha_top = (np.arctan(wingbox_spar()[3][0]-wingbox_spar()[2][0]))/wingbox_spar()[4]
    alpha_bottom = (np.arctan(wingbox_spar()[3][1]-wingbox_spar()[2][1]))/wingbox_spar()[4]
    
    # wingbox = wingbox_NA()
    # wingbox_spar()[0] ---> wingbox['y_NA']
    #Ixx - MOMENT OF INERTIA CONTRIBUTIONS
    front_spar_MOI = t*wingbox_spar()[0]**3 / 12. + wingbox_NA()[2]* (wingbox_NA()[0]-wingbox_NA()[3])**2
    back_spar_MOI = t*wingbox_spar()[1]**3 / 12. + wingbox_NA()[4]* (wingbox_NA()[0]-wingbox_NA()[5])**2
    top_skin_MOI = sum((s1*np.sin(alpha_top)+(np.min(wingbox_spar()[2][0],wingbox_spar()[3][0]) - wingbox_NA()[0]))**2*t*ds)
    bottom_skin_MOI = sum((s1*np.sin(alpha_top))**2*t*ds)
    
    print wingbox_NA()[6], s1
    
    print alpha_top, top_skin_MOI, front_spar_MOI
                                  
    #print back_spar_MOI
    return

wingbox_MOI()

"""




########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#PLOTTING THE AIRFOIL

#plt.figure(figsize = (8.5,6),tight_layout=True)

"""
xfront = 0.15 * np.ones(2)
xback = 0.57 * np.ones(2)
yfront = [-0.04,0.07]
yback = [-0.04,0.09]

plt.plot(xfront,yfront)
plt.plot(xback,yback)


plt.plot(xcoordinates,ycoordinates)


plt.show()
"""






















