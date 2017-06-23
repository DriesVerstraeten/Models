#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:17:57 2017

@author: driesverstraeten
"""
import time
start_time = time.time()

import Init_Parameters as p
import Wing_model as wm
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

t = p.t_skin
c_1 = wm.c_1
c_2 = wm.c_2
c_3 = wm.c_3
c = wm.c
FS_location = 0.15 #front spar location as % of chord
BS_location = 0.575 #front spar location as % of chord

def MOI_section_1():
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    Ixx = []
    Iyy = []
    Ixy = []
    x_span = []
    x_NA = []
    y_NA = []
    LE_area = [] 
    TE_area = [] 
    cross_section_area = []
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    
    for i in range(len(c_1)):
        xcoordinates1 = xcoordinates * c_1[i]
        ycoordinates1 = ycoordinates * c_1[i]
        
        fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
        f1 = np.poly1d(fit1) #upper skin
        f11 = f1.deriv(1)
        
        fit2 = np.polyfit(xcoordinates1[133:],ycoordinates1[133:],5)
        f2 = np.poly1d(fit2) #lower skin
        f22 = f2.deriv(1)

        front_spar = FS_location * c_1[i]
        back_spar = BS_location * c_1[i]

        #sections = 100.             
        #dx = 1./sections             
        x = np.linspace(front_spar,back_spar,len(wm.y))
        dx = (back_spar - front_spar)/len(wm.y)
        area_section = dx * t * 2 #multiplied by two because it is a honeycomb structure
        x_span.append(x)
        '''
        if i == 0:
            plt.plot(xcoordinates1,ycoordinates1)
            plt.axes().set_aspect('equal','datalim')
            plt.plot(x,f1(x),color='k')
            plt.plot(x,f2(x),color='k')
            plt.plot(x,f1(x)-0.008,color='k')
            plt.plot(x,f2(x)+0.008,color='k')
            plt.plot((front_spar,front_spar),(f1(front_spar),f2(front_spar)),color = 'k')
            plt.plot((back_spar,back_spar),(f1(back_spar),f2(back_spar)),color = 'k')
            plt.xlabel('Chord [m]')
            plt.ylabel('Thickness [m]')
            plt.show()
         '''   
        x_LE = np.linspace(0,front_spar,len(wm.y))
        x_TE = np.linspace(back_spar,c_1[i],len(wm.y))

        arc_length_US = integrate.quad(lambda x: np.sqrt(1+f11(x)**2), front_spar, back_spar) #upper skin
        arc_length_LS = integrate.quad(lambda x: np.sqrt(1+f22(x)**2), front_spar, back_spar) #lower skin
        LE_area_sec = wm.dy*(integrate.quad(lambda x_LE: np.sqrt(1+f11(x_LE)**2), 0, front_spar)[0] + integrate.quad(lambda x_LE: np.sqrt(1+f22(x_LE)**2), 0, front_spar)[0])
        TE_area_sec = wm.dy*(integrate.quad(lambda x_TE: np.sqrt(1+f11(x_TE)**2), back_spar, c_1[i])[0] + integrate.quad(lambda x_TE: np.sqrt(1+f22(x_TE)**2), back_spar, c_1[i])[0])                                                       
                                      
        front_spar_length = np.abs(f1(front_spar)) + np.abs(f2(front_spar))
        back_spar_length = np.abs(f1(back_spar)) + np.abs(f2(back_spar))
        
        front_spar_area = (front_spar_length)*t
        back_spar_area = (back_spar_length)*t

        total_area = (arc_length_US[0] + arc_length_LS[0])*2*t + front_spar_area + back_spar_area

        A_crosssection = sum((np.abs(f1(x)) + np.abs(f2(x)))*dx)
        cross_section_area.append(A_crosssection)

        y_NA_sec = (area_section * (sum(f1(x)) + sum(f2(x))) + front_spar_area*(f1(front_spar) + f2(front_spar))/2 + back_spar_area*(f1(back_spar) + f2(back_spar))/2)/ total_area 
        x_NA_sec = (area_section * sum(x) * 2 + front_spar_area * front_spar + back_spar_area * back_spar) / total_area
                   
        x_NA.append(x_NA_sec)
        y_NA.append(y_NA_sec)
        
        Ixx_section = area_section*(sum(f1(x)**2) + sum(f2(x)**2)) + t/12 * (front_spar_length**3 + front_spar_area*(f1(front_spar) + f2(front_spar - y_NA[i]))**2 + back_spar_length**3 + back_spar_area*(f1(back_spar) + f2(back_spar - y_NA[i]))**2)
        Iyy_section = front_spar_area * (front_spar-x_NA[i])**2 + back_spar_area * (back_spar-x_NA[i])**2 + sum(2*area_section*(x-x_NA[i])**2)    #(sum(x**2)-2*x_NA[i]*sum(x)+len(x)*x_NA[i]**2)
        Ixy_section = area_section*(sum((x-x_NA[i])*(f1(x)-y_NA[i])) + sum((x-x_NA[i])*(f2(x)-y_NA[i]))) + front_spar_area * (front_spar - x_NA[i]) * ((f1(front_spar) + f2(front_spar))/2 - y_NA[i]) + back_spar_area * (back_spar - x_NA[i]) * ((f1(back_spar) + f2(back_spar))/2 - y_NA[i])

        Ixx.append(Ixx_section)
        Iyy.append(Iyy_section)
        Ixy.append(Ixy_section)
        LE_area.append(LE_area_sec)
        TE_area.append(TE_area_sec)
        
    return np.array(Ixx), np.array(Iyy), np.array(Ixy), f1, f2, y_NA, x_NA, x, x_span, LE_area, TE_area, cross_section_area


def MOI_section_2():
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    Ixx = []
    Iyy = []
    Ixy = []
    x_span = []
    x_NA = []
    y_NA = []
    LE_area = [] 
    TE_area = [] 
    cross_section_area = []
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    
    for i in range(len(c_2)):
        xcoordinates1 = xcoordinates * c_2[i]
        ycoordinates1 = ycoordinates * c_2[i]
        
        fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
        f1 = np.poly1d(fit1) #upper skin
        f11 = f1.deriv(1)
        
        fit2 = np.polyfit(xcoordinates1[133:],ycoordinates1[133:],5)
        f2 = np.poly1d(fit2) #lower skin
        f22 = f2.deriv(1)
        
        front_spar = FS_location * c_2[i]
        back_spar = BS_location * c_2[i]

        #sections = 100.             
        #dx = 1./sections             
        x = np.linspace(front_spar,back_spar,len(wm.y))
        dx = (back_spar - front_spar)/len(wm.y)
        area_section = dx * t * 2.
        x_span.append(x)

        
        x_LE = np.linspace(0,front_spar,len(wm.y))
        x_TE = np.linspace(back_spar,c_2[i],len(wm.y))

        arc_length_US = integrate.quad(lambda x: np.sqrt(1+f11(x)**2), front_spar, back_spar) #upper skin
        arc_length_LS = integrate.quad(lambda x: np.sqrt(1+f22(x)**2), front_spar, back_spar) #lower skin
        LE_area_sec = wm.dy*(integrate.quad(lambda x_LE: np.sqrt(1+f11(x_LE)**2), 0, front_spar)[0] + integrate.quad(lambda x_LE: np.sqrt(1+f22(x_LE)**2), 0, front_spar)[0])
        TE_area_sec = wm.dy*(integrate.quad(lambda x_TE: np.sqrt(1+f11(x_TE)**2), back_spar, c_2[i])[0] + integrate.quad(lambda x_TE: np.sqrt(1+f22(x_TE)**2), back_spar, c_2[i])[0])
                                      
        front_spar_length = np.abs(f1(front_spar)) + np.abs(f2(front_spar))
        back_spar_length = np.abs(f1(back_spar)) + np.abs(f2(back_spar))
        
        front_spar_area = (front_spar_length)*t
        back_spar_area = (back_spar_length)*t

        total_area = (arc_length_US[0] + arc_length_LS[0])*2*t + front_spar_area + back_spar_area

        A_crosssection = sum((np.abs(f1(x)) + np.abs(f2(x)))*dx)
        cross_section_area.append(A_crosssection)

        y_NA_sec = (area_section * (sum(f1(x)) + sum(f2(x))) + front_spar_area*(f1(front_spar) + f2(front_spar))/2 + back_spar_area*(f1(back_spar) + f2(back_spar))/2)/ total_area 
        x_NA_sec = (area_section * sum(x) * 2 + front_spar_area * front_spar + back_spar_area * back_spar) / total_area
                   
        x_NA.append(x_NA_sec)
        y_NA.append(y_NA_sec)
        
        Ixx_section = area_section*(sum(f1(x)**2) + sum(f2(x)**2)) + t/12 * (front_spar_length**3 + front_spar_area*(f1(front_spar) + f2(front_spar - y_NA[i]))**2 + back_spar_length**3 + back_spar_area*(f1(back_spar) + f2(back_spar - y_NA[i]))**2)
        Iyy_section = front_spar_area * (front_spar-x_NA[i])**2 + back_spar_area * (back_spar-x_NA[i])**2 + sum(2*area_section*(x-x_NA[i])**2)    #(sum(x**2)-2*x_NA[i]*sum(x)+len(x)*x_NA[i]**2)
        Ixy_section = area_section*(sum((x-x_NA[i])*(f1(x)-y_NA[i])) + sum((x-x_NA[i])*(f2(x)-y_NA[i]))) + front_spar_area * (front_spar - x_NA[i]) * ((f1(front_spar) + f2(front_spar))/2 - y_NA[i]) + back_spar_area * (back_spar - x_NA[i]) * ((f1(back_spar) + f2(back_spar))/2 - y_NA[i])

        Ixx.append(Ixx_section)
        Iyy.append(Iyy_section)
        Ixy.append(Ixy_section)
        LE_area.append(LE_area_sec)
        TE_area.append(TE_area_sec)

        
    return np.array(Ixx), np.array(Iyy), np.array(Ixy), f1, f2, y_NA, x_NA, x, x_span, LE_area, TE_area, cross_section_area


def MOI_section_3():
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    Ixx = []
    Iyy = []
    Ixy = []
    x_span = []
    x_NA = []
    y_NA = []
    LE_area = [] 
    TE_area = [] 
    cross_section_area = []
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    
    for i in range(len(c_3)):
        xcoordinates1 = xcoordinates * c_3[i]
        ycoordinates1 = ycoordinates * c_3[i]
        
        fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
        f1 = np.poly1d(fit1) #upper skin
        f11 = f1.deriv(1)
        
        fit2 = np.polyfit(xcoordinates1[133:],ycoordinates1[133:],5)
        f2 = np.poly1d(fit2) #lower skin
        f22 = f2.deriv(1)
        
        #plt.plot(xcoordinates1,ycoordinates1)
        #plt.axes().set_aspect('equal','datalim')
        #plt.show()
        
        front_spar = FS_location * c_3[i]
        back_spar = BS_location * c_3[i]

        #sections = 100.             
        #dx = 1./sections             
        x = np.linspace(front_spar,back_spar,len(wm.y))
        dx = (back_spar - front_spar)/len(wm.y)
        area_section = dx * t * 2.
        x_span.append(x)

        x_LE = np.linspace(0,front_spar,len(wm.y))
        x_TE = np.linspace(back_spar,c_3[i],len(wm.y))

        arc_length_US = integrate.quad(lambda x: np.sqrt(1+f11(x)**2), front_spar, back_spar) #upper skin
        arc_length_LS = integrate.quad(lambda x: np.sqrt(1+f22(x)**2), front_spar, back_spar) #lower skin
        LE_area_sec = wm.dy*(integrate.quad(lambda x_LE: np.sqrt(1+f11(x_LE)**2), 0, front_spar)[0] + integrate.quad(lambda x_LE: np.sqrt(1+f22(x_LE)**2), 0, front_spar)[0])
        TE_area_sec = wm.dy*(integrate.quad(lambda x_TE: np.sqrt(1+f11(x_TE)**2), back_spar, c_3[i])[0] + integrate.quad(lambda x_TE: np.sqrt(1+f22(x_TE)**2), back_spar, c_3[i])[0])
                                      
        front_spar_length = np.abs(f1(front_spar)) + np.abs(f2(front_spar))
        back_spar_length = np.abs(f1(back_spar)) + np.abs(f2(back_spar))
        
        front_spar_area = (front_spar_length)*t
        back_spar_area = (back_spar_length)*t

        total_area = (arc_length_US[0] + arc_length_LS[0])*2*t + front_spar_area + back_spar_area

        A_crosssection = sum((np.abs(f1(x)) + np.abs(f2(x)))*dx)
        cross_section_area.append(A_crosssection)

        y_NA_sec = (area_section * (sum(f1(x)) + sum(f2(x))) + front_spar_area*(f1(front_spar) + f2(front_spar))/2 + back_spar_area*(f1(back_spar) + f2(back_spar))/2)/ total_area 
        x_NA_sec = (area_section * sum(x) * 2 + front_spar_area * front_spar + back_spar_area * back_spar) / total_area
                   
        x_NA.append(x_NA_sec)
        y_NA.append(y_NA_sec)

        Ixx_section = area_section*(sum(f1(x)**2) + sum(f2(x)**2)) + t/12 * (front_spar_length**3 + front_spar_area*(f1(front_spar) + f2(front_spar - y_NA[i]))**2 + back_spar_length**3 + back_spar_area*(f1(back_spar) + f2(back_spar - y_NA[i]))**2)
        Iyy_section = front_spar_area * (front_spar-x_NA[i])**2 + back_spar_area * (back_spar-x_NA[i])**2 + sum(2*area_section*(x-x_NA[i])**2)    #(sum(x**2)-2*x_NA[i]*sum(x)+len(x)*x_NA[i]**2)
        Ixy_section = area_section*(sum((x-x_NA[i])*(f1(x)-y_NA[i])) + sum((x-x_NA[i])*(f2(x)-y_NA[i]))) + front_spar_area * (front_spar - x_NA[i]) * ((f1(front_spar) + f2(front_spar))/2 - y_NA[i]) + back_spar_area * (back_spar - x_NA[i]) * ((f1(back_spar) + f2(back_spar))/2 - y_NA[i])

        Ixx.append(Ixx_section)
        Iyy.append(Iyy_section)
        Ixy.append(Ixy_section)
        LE_area.append(LE_area_sec)
        TE_area.append(TE_area_sec)
        
    return np.array(Ixx), np.array(Iyy), np.array(Ixy), f1, f2, y_NA, x_NA, x, x_span, LE_area, TE_area, cross_section_area

Ixx_1, Iyy_1, Ixy_1, f1_1, f2_1, y_NA_1, x_NA_1, x_1, x_span_1, LE_area_1, TE_area_1, ribs_area_1 = MOI_section_1()
Ixx_2, Iyy_2, Ixy_2, f1_2, f2_2, y_NA_2, x_NA_2, x_2, x_span_2, LE_area_2, TE_area_2, ribs_area_2 = MOI_section_2()
Ixx_3, Iyy_3, Ixy_3, f1_3, f2_3, y_NA_3, x_NA_3, x_3, x_span_3, LE_area_3, TE_area_3, ribs_area_3 = MOI_section_3()

def wingbox_MOI_total():
    Ixx = np.hstack((Ixx_1, Ixx_2, Ixx_3)) 
    Iyy = np.hstack((Iyy_1, Iyy_2, Iyy_3))
    Ixy = np.hstack((Ixy_1, Ixy_2, Ixy_3))
    x_NA = np.hstack((x_NA_1, x_NA_2, x_NA_3))
    y_NA = np.hstack((y_NA_1, y_NA_2, y_NA_3))
    x_span = np.vstack((x_span_1, x_span_2, x_span_3))
    LE_area = sum(np.hstack((LE_area_1, LE_area_2, LE_area_3)))
    TE_area = sum(np.hstack((TE_area_1, TE_area_2, TE_area_3)))
    ribs_area = np.hstack((ribs_area_1[0], ribs_area_2[0], ribs_area_3[0]))
    total_volume_wings = 2*wm.dy*(sum(np.hstack((ribs_area_1, ribs_area_2, ribs_area_3))))
    
    return Ixx, Iyy, Ixy, x_NA, y_NA, x_span, LE_area, TE_area, ribs_area, total_volume_wings


Ixx, Iyy, Ixy, x_NA, y_NA, x_span, LE_area, TE_area, ribs_area, total_volume_wings = wingbox_MOI_total()

print "MOI:", total_volume_wings, "m^3"


print("MOI --- %s seconds ---" % (time.time() - start_time))

#plt.close()
#plt.plot(wm.y,Ixx, color = 'r')
#plt.plot(wm.y,Iyy, color = 'b')
#plt.plot(wm.y,Ixy, color = 'g')
#plt.plot(wm.y,Ixx*Iyy - Ixy**2, 'r')
#plt.show()


'''



########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
#PLOTTING THE AIRFOIL

#plt.figure(figsize = (8.5,6),tight_layout=True)


xfront = 0.15 * np.ones(2)
xback = 0.57 * np.ones(2)
yfront = [-0.04,0.07]
yback = [-0.04,0.09]

plt.plot(xfront,yfront)
plt.plot(xback,yback)





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

'''



