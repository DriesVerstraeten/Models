#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:56:43 2017

@author: driesverstraeten
"""
import MOI as mi
import Wing_model as wm
import Init_Parameters as p
import numpy as np
import matplotlib.pyplot as plt
import time 
start_time = time.time()


c_1 = wm.c_1
c_2 = wm.c_2
c_3 = wm.c_3
Mx = wm.wing_moment(wm.CL,p.rho_0,p.V_cruise)[1]
Ixx, Iyy, Ixy, x_NA, y_NA, x_span, LE_area, TE_area, ribs_area, total_volume_wings = mi.wingbox_MOI_total()

def piecewise_poly_1():
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    y_position_US = []
    y_position_LS = []
    y_position_FS_u = []
    y_position_FS_l = []
    y_position_BS_u = []
    y_position_BS_l = []
    
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    for i in range(len(c_1)):
        xcoordinates1 = xcoordinates * c_1[i]
        ycoordinates1 = ycoordinates * c_1[i]
        
        fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
        f1 = np.poly1d(fit1)
        y_position_US.append(f1(x_span[i]))
        y_position_FS_u.append(f1(x_span[i][0]))
        y_position_BS_u.append(f1(x_span[i][-1]))
        
        fit2 = np.polyfit(xcoordinates1[133:],ycoordinates1[133:],5)
        f2 = np.poly1d(fit2)
        y_position_LS.append(f2(x_span[i]))
        y_position_FS_l.append(f2(x_span[i][0]))
        y_position_BS_l.append(f2(x_span[i][-1]))
    
    return np.array(y_position_US), np.array(y_position_LS), np.array(y_position_FS_u), np.array(y_position_BS_u), np.array(y_position_FS_l), np.array(y_position_BS_l)




def piecewise_poly_2():
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    y_position_US = []
    y_position_LS = []
    y_position_FS_u = []
    y_position_FS_l = []
    y_position_BS_u = []
    y_position_BS_l = []
    
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    for i in range(len(c_2)):
        xcoordinates1 = xcoordinates * c_2[i]
        ycoordinates1 = ycoordinates * c_2[i]
        
        fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
        f1 = np.poly1d(fit1)
        y_position_US.append(f1(x_span[i]))
        y_position_FS_u.append(f1(x_span[i][0]))
        y_position_BS_u.append(f1(x_span[i][-1]))
        
        fit2 = np.polyfit(xcoordinates1[133:],ycoordinates1[133:],5)
        f2 = np.poly1d(fit2)
        y_position_LS.append(f2(x_span[i]))
        y_position_FS_l.append(f2(x_span[i][0]))
        y_position_BS_l.append(f2(x_span[i][-1]))
    
    return np.array(y_position_US), np.array(y_position_LS), np.array(y_position_FS_u), np.array(y_position_BS_u), np.array(y_position_FS_l), np.array(y_position_BS_l)




def piecewise_poly_3():
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    y_position_US = []
    y_position_LS = []
    y_position_FS_u = []
    y_position_FS_l = []
    y_position_BS_u = []
    y_position_BS_l = []
    
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    for i in range(len(c_3)):
        xcoordinates1 = xcoordinates * c_3[i]
        ycoordinates1 = ycoordinates * c_3[i]
        
        fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
        f1 = np.poly1d(fit1)
        y_position_US.append(f1(x_span[i]))
        y_position_FS_u.append(f1(x_span[i][0]))
        y_position_BS_u.append(f1(x_span[i][-1]))
        
        fit2 = np.polyfit(xcoordinates1[133:],ycoordinates1[133:],5)
        f2 = np.poly1d(fit2)
        y_position_LS.append(f2(x_span[i]))
        y_position_FS_l.append(f2(x_span[i][0]))
        y_position_BS_l.append(f2(x_span[i][-1]))
    
    return np.array(y_position_US), np.array(y_position_LS), np.array(y_position_FS_u), np.array(y_position_BS_u), np.array(y_position_FS_l), np.array(y_position_BS_l)

poly1_out = piecewise_poly_1()
poly2_out = piecewise_poly_2()
poly3_out = piecewise_poly_3()
y_position_US = np.vstack((poly1_out[0], poly2_out[0], poly3_out[0]))
y_position_LS = np.vstack((poly1_out[1], poly2_out[1], poly3_out[1]))
y_position_FS_u = np.hstack((poly1_out[2], poly2_out[2], poly3_out[2]))
y_position_BS_u = np.hstack((poly1_out[3], poly2_out[3], poly3_out[3]))
y_position_FS_l = np.hstack((poly1_out[4], poly2_out[4], poly3_out[4]))
y_position_BS_l = np.hstack((poly1_out[5], poly2_out[5], poly3_out[5]))




def wingbox_bending_stress():
    sigma_bending_US = np.ones((len(wm.y),len(Ixx)))
    sigma_bending_LS = np.ones((len(wm.y),len(Ixx)))
    sigma_bending_FS = np.ones((len(wm.y),len(Ixx)))
    sigma_bending_BS = np.ones((len(wm.y),len(Ixx)))
    
    for i in range(len(wm.y)):
        x = x_span[i] - x_NA[i]
        y_US = y_position_US[i] - y_NA[i]
        y_LS = y_position_LS[i] - y_NA[i]
        y_FS = np.linspace(y_position_FS_u[i],y_position_FS_l[i],len(wm.y))-y_NA[i]
        y_BS = np.linspace(y_position_BS_u[i],y_position_BS_l[i],len(wm.y))-y_NA[i]
        
        
        sigma_bending_US[i,:] = - Mx[i] / (Ixx[i] * Iyy[i] - Ixy[i]**2) * (Iyy[i] * y_US - Ixy[i] * x)
        sigma_bending_LS[i,:] = - Mx[i] / (Ixx[i] * Iyy[i] - Ixy[i]**2) * (Iyy[i] * y_LS - Ixy[i] * x)
        sigma_bending_FS[i,:] = - Mx[i] / (Ixx[i] * Iyy[i] - Ixy[i]**2) * (Iyy[i] * y_FS - Ixy[i] * x[0])
        sigma_bending_BS[i,:] = - Mx[i] / (Ixx[i] * Iyy[i] - Ixy[i]**2) * (Iyy[i] * y_BS - Ixy[i] * x[-1])
    
    return sigma_bending_US, sigma_bending_LS, sigma_bending_FS, sigma_bending_BS

sigma_bending_US, sigma_bending_LS, sigma_bending_FS, sigma_bending_BS = wingbox_bending_stress()

maxpoints = np.ones(len(wm.y))
for i in range(len(wm.y)):
    
    maxpoints[i] = np.max(sigma_bending_US[i])

#plt.plot(wm.y, y_US[0])
#plt.plot(x_span[0],sigma_bending_US[0], color='r')
#plt.plot(x_span[0],sigma_bending_LS[0], color='b')
#plt.plot(np.ones(len(wm.y))*x_span[0][0],sigma_bending_FS[0], color='y')
#plt.plot(np.ones(len(wm.y))*x_span[0][-1],sigma_bending_BS[0], color='m')
plt.plot(wm.y,maxpoints, color='b')
plt.show()






'''
c = wm.c
f1 = mi.wingbox_MOI()[3]
f2 = mi.wingbox_MOI()[4]
y_NA = mi.wingbox_MOI()[5]
x_NA = mi.wingbox_MOI()[6]
x = mi.wingbox_MOI()[7]
x_span1 = mi.wingbox_MOI()[8]
x_span = np.ones(len(x_span1)) #x_span is the x locations at every spanwise piece
Mx = wm.wing_moment_9g(wm.CL_9g,p.rho_0,p.V_cruise)[1]
Ixx = mi.wingbox_MOI()[0]
Iyy = mi.wingbox_MOI()[1]
Ixy = mi.wingbox_MOI()[2]

#print Mx[0],Ixx[0],Iyy[0],Ixy[0],x_NA[0],y_NA[0]

def wingbox_extreme_pos():
    longest_US = np.sqrt((f1(x)-y_NA[-1])**2 + (x-x_NA[-1])**2)
    longest_LS = np.sqrt((f2(x)-y_NA[-1])**2 + (x-x_NA[-1])**2)

    longest_pos_US = [i for i,p in enumerate(longest_US) if p == np.max(longest_US)]
    longest_pos_LS = [i for i,p in enumerate(longest_LS) if p == np.max(longest_LS)]
    
    #print longest_pos_US, longest_pos_LS
    #print longest_US[longest_pos_US], longest_LS[longest_pos_LS]
    
    if longest_US[longest_pos_US] > longest_LS[longest_pos_LS]:
        x_location = longest_pos_US[0]
        f = f1
    else: 
        x_location = longest_pos_LS[0]
        f = f2

    return longest_US, longest_LS, x_location, f

wingbox_extreme_pos()
longest_US = wingbox_extreme_pos()[0]
longest_LS = wingbox_extreme_pos()[1]
x_location = wingbox_extreme_pos()[2]
f = wingbox_extreme_pos()[3]
#print f




for i in range(len(x_span)):
    x_span[i] = x_span1[i][x_location]


#print x_span1[0][x_location]



def extreme_positions():
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    y_position_US = []
    y_position_LS = []
    
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    for i in range(len(c)):
        xcoordinates1 = xcoordinates * c[i]
        ycoordinates1 = ycoordinates * c[i]
        
        fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
        f1 = np.poly1d(fit1)
        y_position_US.append(f1(x_span[i]))
        
        #print f1(x_span[x_location])
        
        fit2 = np.polyfit(xcoordinates1[133:256],ycoordinates1[133:256],5)
        f2 = np.poly1d(fit2)
        y_position_LS.append(f2(x_span[i]))
    
    return np.array(y_position_US), np.array(y_position_LS)

y_position_US = extreme_positions()[0]
y_position_LS = extreme_positions()[1]





def wingbox_bending_stress():
    Mx = wm.wing_moment_9g(wm.CL_9g,p.rho_0,p.V_cruise)[1]
    Ixx = mi.wingbox_MOI()[0]
    Iyy = mi.wingbox_MOI()[1]
    Ixy = mi.wingbox_MOI()[2]
    #x_location = wingbox_extreme_pos()[2]
    f = wingbox_extreme_pos()[3]
    
    max_bending_stress = - Mx[1:] / (Ixx*Iyy - Ixy**2) * (Iyy*(y_position_LS-y_NA) - Ixy*((x_span-x_NA))) #y_NA changes with the spanwise sections, so change this!!
    max_bending_stressu = - Mx[1:] / (Ixx*Iyy - Ixy**2) * (Iyy*(y_position_US-y_NA) - Ixy*((x_span-x_NA)))
    #print max_bending_stress, max_bending_stressu
    
    plt.plot(wm.y,x_span-x_NA, ' og' )
    plt.show()
    
    return max_bending_stress, f, max_bending_stressu

wingbox_bending_stress()


plt.plot(wm.y,wingbox_bending_stress()[2])
plt.show()

maxbendstress = wingbox_bending_stress()[0]
maxbendstressu = wingbox_bending_stress()[2]
'''
print("--- %s seconds ---" % (time.time() - start_time))







