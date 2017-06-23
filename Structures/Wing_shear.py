# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 13:57:22 2017

@author: 123
"""
import Init_Parameters as p
import numpy as np
import matplotlib.pyplot as plt
import math
import MOI as moi
#import Wing_model as wing


dy = 0.01 #small spanwise section

r2 = .95
r3 = .8
r4 = .4

def ChordBase(r2,r3,r4,b,S):

    return S/((r4/2.0+r3+r2+.5)*(b/3.0))

def GetChords(r2,r3,r4,b,S):

    C1 = ChordBase(r2,r3,r4,b,S)
    C2 = C1*r2
    C3 = C1*r3
    C4 = C1*r4

    return C1,C2,C3,C4


C1,C2,C3,C4 = GetChords(r2,r3,r4,p.b,p.S)

bsection = p.b/6.

c_1 = np.linspace(C1, C2, bsection*1./dy)
c_2 = np.linspace(C2, C3, bsection*1./dy)
c_3 = np.linspace(C3, C4, bsection*1./dy)
c = np.hstack((c_1,c_2,c_3))

y = np.linspace(0,p.b/2.+dy, len(c)) #spanwise location of section
y1 = np.linspace(0,p.b/2.+2*dy, len(c))

Ixx = moi.wingbox_MOI_total()[0] 
Iyy = moi.wingbox_MOI_total()[1]
Ixy = moi.wingbox_MOI_total()[2]
y_NA = moi.wingbox_MOI_total()[4]
x_NA = moi.wingbox_MOI_total()[3]

def shear_box(dx,i,t,fh,T): #shear stress calc. dx dy are grid numbers (use 500 for each or so). start end are %cord where wingbox is (0.1 and 0.575 for wing). i is the location where you calculate shear (so if its 499 for dx=500 its at the root. or tip if you defined it like that). t thickness, rho is rho, fh is the force on the wing
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    x_x = np.linspace(0.15*c[i],0.575*c[i],dx)
    x_x_inv = np.linspace(0.575*c[i],0.15*c[i],dx)
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
    for j in range(len(airfoil_coordinates)):
        xcoordinates[j] = airfoil_coordinates[j][0]
        ycoordinates[j] = airfoil_coordinates[j][1] 
    xcoordinates1 = xcoordinates * c[i]
    ycoordinates1 = ycoordinates * c[i]
    fit1 = np.polyfit(xcoordinates1[0:133],ycoordinates1[0:133],5)
    f1 = np.poly1d(fit1) 
    f11 = f1.deriv(1)
    fit2 = np.polyfit(xcoordinates1[133:256],ycoordinates1[133:256],5)
    f2 = np.poly1d(fit2)
    f22 = f2.deriv(1)
    x_y_up = np.zeros(dx)
    x_y_down = np.zeros(dx)
    for j in range (0,dx):
        x_y_up[j] = f1(x_x[j])
        x_y_down[j] = f2(x_x[j])
    x_y_right = np.linspace(x_y_up[dx-1],x_y_down[dx-1],dx) 
    x_y_left = np.linspace(x_y_down[0],x_y_up[0],dx) 
    slope1 = np.zeros(dx)
    slope2 = np.zeros(dx)
    for j in range (0,dx):
        slope1[j] = ((-0.15*c[i]+0.575*c[i])/dx)/(math.cos(math.atan(f11(x_x[j]))))
        slope2[j] = ((-0.15*c[i]+0.575*c[i])/dx)/(math.cos(math.atan(f22(x_x[j]))))
    q1 = np.zeros(dx)
    q2 = np.zeros(dx)
    q3 = np.zeros(dx)
    q4 = np.zeros(dx)
    q1[x_NA[i]] = 0
    for j in range (int(x_NA[i+1]*dx/c[i]),dx):
        q1[j] = q1[j-1] - (fh*slope1[x_x[j]]*t*((x_y_up[j]-y_NA[i])*Iyy[i]-(x_x[j]-x_NA[i])*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    q2[0] = q1[dx-1]
    for j in range (1,dx):
        q2[j] = q2[j-1] - (fh*(abs(x_y_down[dx-1]-x_y_up[dx-1])/dx)*t*((x_y_right[j]-y_NA[i])*Iyy[i]-(x_x[dx-1]-x_NA[i])*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    q3[0] = q2[dx-1]
    for j in range (1,dx):
        q3[j] = q3[j-1] - (fh*slope2[x_x[dx-j]]*t*((x_y_down[dx-j]-y_NA[i])*Iyy[i]-(x_x[dx-j]-x_NA[i])*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    q4[0] = q3[dx-1]
    for j in range (1,dx):
        q4[j] = q4[j-1] - (fh*(abs(x_y_down[0]-x_y_up[0])/dx)*t*((x_y_left[j]-y_NA[i])*Iyy[i]-(x_x[0]-x_NA[i])*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    q1[0] = q4[dx-1]
    for j in range (1,int(x_NA[i]*dx/c[i])):
        q1[j] = q1[j-1] - (fh*slope1[x_x[j]]*t*((x_y_up[j]-y_NA[i])*Iyy[i]-(x_x[j]-x_NA[i])*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    encl_area = 0
    for j in range (0,dx):
        encl_area = encl_area + (x_y_up[j]-x_y_down[j])/dx
    qt = 2 * encl_area * T
    q1 = q1 + qt
    q2 = q2 + qt
    q3 = q3 + qt
    q4 = q4 + qt
    mx = max(abs(max(q1)),abs(min(q1)),abs(max(q2)),abs(min(q2)),abs(max(q3)),abs(min(q3)),abs(max(q4)),abs(min(q4)))
    plt.figure(figsize=(19,5), tight_layout = True)
    plt.suptitle('Shear flow in the wingbox')
    plt.subplot(141)
    plt.plot(x_x,q1)
    plt.ylabel('Shear flow top skin, N/m^2')
    plt.xlabel('x - Location')
    plt.subplot(142)
    plt.plot(x_y_right,q2)
    plt.ylabel('Shear flow right spar, N/m^2')
    plt.xlabel('y - Location')
    plt.subplot(143)
    plt.plot(x_x_inv,q3)
    plt.ylabel('Shear flow bottom skin, N/m^2')
    plt.xlabel('x - Location')
    plt.subplot(144)
    plt.plot(x_y_right,q4)
    plt.ylabel('Shear flow left spar, N/m^2')
    plt.xlabel('y - Location')
    plt.show()
    return mx
 

def highest_shear_box(dx,t,fh,rib1,rib2,T):
    y = np.linspace(0,p.b/2,300)
    stress = np.zeros(300)
    for i in range (0,300):
        stress[i] = shear_box(dx,i,t,fh,rib1,rib2,T)/t
    plt.figure(figsize=(7,5))
    plt.suptitle('Highest shear stress in the tailbox')
    plt.plot(y,stress)
    plt.ylabel('Shear stress, N/m^2')
    plt.xlabel('x - Location')
    plt.show()
    a = max(stress)
    b = min(stress)
    if abs(a) >= abs(b):
        c = a
    else:
        c = b
    return c