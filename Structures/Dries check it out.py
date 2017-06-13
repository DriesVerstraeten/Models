# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 13:57:22 2017

@author: 123
"""

def shear_box(dx,dy,start,end,i,t,rho,fh): #shear stress calc. dx dy are grid numbers (use 500 for each or so). start end are %cord where wingbox is (0.1 and 0.575 for wing). i is the location where you calculate shear (so if its 499 for dx=500 its at the root. or tip if you defined it like that). t thickness, rho is rho, fh is the force on the wing
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    x_x = np.linspace(start*cord at i,i),end*cord at i,dx)#insert formula for cord at i
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
    for j in range(len(airfoil_coordinates)):
        xcoordinates[j] = airfoil_coordinates[j][0]
        ycoordinates[j] = airfoil_coordinates[j][1] 
    xcoordinates1 = xcoordinates * cord at i#insert formula or function for it here
    ycoordinates1 = ycoordinates * cord at i
    fit1 = np.polyfit(xcoordinates1[0:66],ycoordinates1[0:66],5)#change number back!
    f1 = np.poly1d(fit1)
    f11 = f1.deriv(1)
    fit2 = np.polyfit(xcoordinates1[68:133],ycoordinates1[68:133],5)
    f2 = np.poly1d(fit2)
    f22 = f2.deriv(1)
    x_y_up = np.zeros(dx)
    x_y_down = np.zeros(dx)
    for j in range (0,dx):
        x_y_up[j] = f1(x_x[j])
        x_y_down[j] = f2(x_x[j])
    x_y_right = np.linspace(x_y_up[dx-1],x_y_down[dx-1],dx) 
    x_y_left = np.linspace(x_y_down[0],x_y_up[0],dx) 
    Ixx = wingbox_MOI(dx,start,end,t,rho)[0] #check for inputs, i have modified function for myself a bit, yours might be different
    Iyy = wingbox_MOI(dx,start,end,t,rho)[1]
    Ixy = wingbox_MOI(dx,start,end,t,rho)[2]
    y_NA = wingbox_MOI(dx,start,end,t,rho)[8]
    x_NA = wingbox_MOI(dx,start,end,t,rho)[9]
    slope1 = np.zeros(dx)
    slope2 = np.zeros(dx)
    for j in range (0,dx):
        slope1[j] = ((-start*c_tail(dy,i)+end*c_tail(dy,i))/dx)/(math.cos(math.atan(f11(x_x[j]))))#replace c_tail by cord values
        slope2[j] = ((-start*c_tail(dy,i)+end*c_tail(dy,i))/dx)/(math.cos(math.atan(f22(x_x[j]))))
    q1 = np.zeros(dx)
    q2 = np.zeros(dx)
    q3 = np.zeros(dx)
    q1[0] = 0
    for j in range (1,dx):
        q1[j] = q1[j-1] - (fh*slope1[j]*t*((x_y_up[j]-y_NA)*Iyy[i]-(x_x[j]-x_NA)*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    q2[0] = q1[dx-1]
    for j in range (1,dx):
        q2[j] = q2[j-1] - (fh*(abs(x_y_down[dx-1]-x_y_up[dx-1])/dx)*t*((x_y_right[j]-y_NA)*Iyy[i]-(x_x[dx-1]-x_NA[i])*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    q3[dx-1] = q2[dx-1]
    for j in range (dx-2,0):
        q3[j] = q3[j+1] - (fh*slope2[j]*t*((x_y_down[j]-y_NA)*Iyy[i]-(x_x[j]-x_NA)*Ixy[i]))/(Ixx[i]*Iyy[i]-Ixy[i]**2)
    for j in range (1,dx): in case of additional rip uncoment this
        q[j+3*dx] = q[3*dx+j-1] + fh*(abs(x_y_down[0]-x_y_up[0])/dx)*t*(x_y_down[0]+j*(abs(x_y_down[0]-x_y_up[0])/dx)/Ixx[i]
    plt.figure(figsize=(19,5))
    plt.suptitle('Shear flow in the tailbox')
    plt.subplot(131)
    plt.plot(x_x,q1)
    plt.ylabel('Shear flow top skin, N/m^2')
    plt.xlabel('x - Location')
    plt.subplot(132)
    plt.plot(x_y_right,q2)
    plt.ylabel('Shear flow right spar, N/m^2')
    plt.xlabel('x - Location')
    plt.subplot(133)
    plt.plot(x_x,q3)
    plt.ylabel('Shear flow bottom skin, N/m^2')
    plt.xlabel('x - Location')
    plt.show()
    return max(q2)

def highest_shear_box(dx,dy,start,end,t,rho,fh):
    y = np.linspace(0,p.b_ht/2,dy)
    stress = np.zeros(dy)
    for i in range (0,dy):
        stress[i] = shear_box(dx,dy,start,end,i,t,rho,fh)/t
    plt.figure(figsize=(7,5))
    plt.suptitle('Highest shear stress in the tailbox')
    plt.plot(y,stress)
    plt.ylabel('Bending stress, N/m^2')
    plt.xlabel('x - Location')
    plt.show()
    a = max(stress)
    b = min(stress)
    if abs(a) >= abs(b):
        c = a
    else:
        c = b
    return c