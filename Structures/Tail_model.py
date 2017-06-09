import Init_Parameters as p
import Material_properties as mat
import numpy as np
import matplotlib.pyplot as plt
import math
"""
INPUTS OF DIMENSIONS ARE HERE FOR NOW. everything in METERS
"""
#WINGBOX DIMENSIONS
box_t = 0.0005 #m, thickness. same everywhere for now
box_b = 0.25*p.cr_ht #m, width at rootapply scaling for other points. scales linearly with taper
box_h = 0.1*p.cr_ht #m, height at root. apply scaling for other points. scales linearly with taper

#boom design
b_rr = 0.9 #root radius
b_rt = 0.4 #tail radius
b_t = 0.001 #thickness
b_l = 2 #length

#Vertical tail dimensions
v_r = p.cr_ht#root cord length
v_t = p.ct_ht#tip cord length
v_b = p.b_ht / 2#vertical tail span
v_w = p.W_ht#weight of the vertical tail. gonna be recalculated a bit later

"""
MAIN FUNCTIONS ARE HERE. inputs are yada yada

TO-DO list: integrate MOI from Dries file, integrate shape of the wingbox for shear distribution with Dries help
Output min/max stress and say if the chosen material can withstand the load. Input material number instead of just rho.
Make master-function.
"""


"""
Stuff

"""
def y_tail(dy,i): #from tip to root, verified
    y_tail = np.linspace(0,p.b_ht/2,dy)
    return y_tail[i]

def y_tail_inv(dy,i): #from root to tip, verified
    y_tail_inv = np.linspace(p.b_ht/2,0,dy)
    return y_tail_inv[i]

def c_tail(dy,i): #from tip to root, works
    c_tail = np.linspace(p.ct_ht,p.cr_ht,dy+1)
    return c_tail[i]

def tail_force(dy,f,i): #from tip to root, works
    w_distr = f * c_tail(dy,i) / c_tail(dy,dy/2)/dy
    return w_distr

def slope_c(dy): #slope of cord, works
    slope = (p.ct_ht+p.cr_ht)/dy
    return slope

def tail_distr(dy,f,i): #loading on a piece of tail (per meter)
    distr = f*c_tail(dy,i) / c_tail(dy,dy/2)*dy/(p.b_ht/2)
    return distr 

def tail_shear(dy,f,j): #works, shear at each point
    shear = np.zeros(j+1)
    shear[0] = tail_force(dy,f,0)
    for i in range (1,j+1):
        shear[i] = shear[i-1] + tail_force(dy,f,i)
    return shear[j]

def tail_moment(dy,f,j): #works
    moment = 0
    for i in range (0,j+1):
        moment = tail_force(dy,f,i) * (y_tail(dy,j)-y_tail(dy,i)) + moment
    return moment
        
def tail_plots(dy,f): #works
    distr = np.zeros(dy)
    y = np.zeros(dy)
    shear = np.zeros(dy)
    moment = np.zeros(dy)
    for i in range (0,dy):
        distr[i] = tail_distr(dy,f,i)
    for i in range (0,dy):
        y[i] = y_tail(dy,i)
    for i in range (0,dy):
        shear[i] = tail_shear(dy,f,i)
    for i in range (0,dy):
        moment[i] = tail_moment(dy,f,i)
    plt.figure(figsize=(19,5))
    plt.suptitle('Load-shear-moment diagrams for horisontal tail')
    plt.subplot(131)
    plt.plot(y, distr)
    plt.ylabel('Horizontal tail load, N')
    plt.xlabel('Location from tip to root, m')
    plt.subplot(132)
    plt.plot(y, shear)
    plt.ylabel('Horizontal tail shear, N')
    plt.xlabel('Location from tip to root, m')
    plt.subplot(133)
    plt.plot(y, moment)
    plt.ylabel('Horizontal tail moment, Nm')
    plt.xlabel('Location from tip to root, m')
    plt.show()
    return

def tail_shear_stress(dy,dx,dz,i,f,Ixx,Ixy,Iyy,t):#looks like it works
    x_12 = np.linspace(0, (1-y_tail_inv(dy,i) * slope_c(dy)) * box_b/2,dx)
    z_23 = np.linspace(((1 - y_tail_inv(dy,i) * slope_c(dy)) * box_h)/2, -1*(1-y_tail_inv(dy,i)*slope_c(dy))*0.1*p.cr_ht/2, dz)
    x_34 = np.linspace((1 - y_tail_inv(dy,i) * slope_c(dy)) * box_b/2,0,dx)
    q_12 = np.zeros(dx)
    q_23 = np.zeros(dy)
    q_34 = np.zeros(dx)
    for j in range (0,dx):
        q_12[j] = -tail_shear(dy,f,i) * t * z_23[0] * x_12[j] / Ixx[i]
    for j in range (0,dz):
        q_23[j] = q_12[dx-1] + box_t * tail_shear(dy,f,i) * x_12[dx-1] * ((z_23[0])**2 - z_23[j]**2)/(4*Ixx[i])
    for j in range (0,dx):
        q_34[j] = q_23[dz-1] - tail_shear(dy,f,i) * box_t * z_23[dz-1] * x_12[j] / Ixx[i]
    tau_12 = np.zeros(dx)
    tau_23 = np.zeros(dx)
    tau_34 = np.zeros(dx)
    for j in range (0,dx):
        tau_12[j] = q_12[j] / t
    for j in range (0,dx):
        tau_23[j] = q_23[j] / t
    for j in range (0,dx):
        tau_34[j] = q_34[j] / t
    plt.figure(figsize=(19,10))
    plt.suptitle('Shear distribution in horisontal tail box')
    plt.subplot(231)
    plt.plot(x_12, q_12)
    plt.ylabel('Shear flow, N/m')
    plt.xlabel('Location (top plate from middle to right edge), m')
    plt.subplot(232)
    plt.plot(q_23, z_23)
    plt.ylabel('Location (right plate from bottom to top), m')
    plt.xlabel('Shear flow, N/m')
    plt.subplot(233)
    plt.plot(x_34, q_34)
    plt.ylabel('Shear flow, N/m')
    plt.xlabel('Location (bottom plate from middle to right), m')
    plt.subplot(234)
    plt.plot(x_12, tau_12)
    plt.ylabel('Shear stress, N/m^2')
    plt.xlabel('Location (top plate from middle to right edge), m')
    plt.subplot(235)
    plt.plot(z_23, tau_23)
    plt.ylabel('Location (right plate from bottom to top), m')
    plt.xlabel('Shear stress, N/m^2')
    plt.subplot(236)
    plt.plot(x_34, tau_34)
    plt.ylabel('Shear stress, N/m^2')
    plt.xlabel('Location (bottom plate from middle to right), m')
    plt.show()   
    return

def x_boom(dx,i): #works
    x = np.linspace(0,b_l,dx)
    return x[i]

def r_boom(dx,i): # from tip to root, works
    slope = (b_rr-b_rt)/dx
    r = b_rt + slope*i
    return r

def W_boom(rho):#works
    w = rho * (math.pi * (b_rr+b_rt)**2/4 * b_l - math.pi * (b_rr+b_rt-2*b_t)**2/4*b_l)
    return w

def w_boom(dx,i):# tip to root, works
    mid_w = W_boom(mat.rho[0]) / dx
    w_distr = mid_w * r_boom(dx,i) / r_boom(dx,dx/2)
    return w_distr

def force_boom(dx,g,i): #works
    force = g*p.g*w_boom(dx,i)
    return force

def distr_boom(dx,g,i,f):#works
    distr = (g*p.g*w_boom(dx,i)*dx + f)/(b_l)
    return distr

def shear_boom(dx,g,j,f):#works
    distr = np.zeros(j+1)
    for i in range (1,j+1):
        distr[i] = force_boom(dx,g,i) 
    shear_j = force_boom(dx,g,0)+ f + np.sum(distr)
    return shear_j

def moment_boom(dx,f,g,j):#works
    moment = 0
    for i in range (0,j+1):
        moment = force_boom(dx,g,i) * (x_boom(dx,j)-x_boom(dx,i)) + moment
    moment = f*x_boom(dx,j) + moment
    return moment

def boom_plots(dx,f,g):#works
    distr = np.zeros(dx)
    x = np.zeros(dx)
    shear = np.zeros(dx)
    moment = np.zeros(dx)
    for i in range (0,dx):
        distr[i] = distr_boom(dx,g,i,f)
    for i in range (0,dx):
        x[i] = x_boom(dx,i)
    shear[0] = shear_boom(dx,g,0,f)
    for i in range (1,dx):
        shear[i] = shear_boom(dx,g,i,f)
    for i in range (0,dx):
        moment[i] = moment_boom(dx,f,g,i)
    plt.figure(figsize=(19,5))
    plt.suptitle('Load-shear-moment diagram for the mounting boom by horisontal tail')
    plt.subplot(131)
    plt.plot(x, distr)
    plt.ylabel('Load, N')
    plt.xlabel('Location, m')
    plt.subplot(132)
    plt.plot(x, shear)
    plt.ylabel('Shear, N')
    plt.xlabel('Location, m')
    plt.subplot(133)
    plt.plot(x, moment)
    plt.ylabel('Moment, Nm')
    plt.xlabel('Location, m')
    plt.show()
    return    
    
def force_vert(dz,f):#works not anymore ARRRRRRRRRRRRRRRGH
    v_c = np.linspace(v_r,v_t,dz)
    mid_w = (f+0.00) / dz
    w_distr = mid_w * v_c / v_c[dz/2]
    return np.sum(w_distr)

def boom_plots_vertical(dx,f):#works
    x = np.zeros(dx)
    for i in range (0,dx):
        x[i] = x_boom(dx,i)
    moment = np.zeros(dx)
    shear = force_vert(1000,f)
    moment = shear * x
    plt.figure(figsize=(9,5))
    plt.suptitle('Moment diagram for the mounting boom by vertical tail')
    plt.plot(x, moment)
    plt.ylabel('Moment, Nm')
    plt.xlabel('Location, m')
    plt.show()
    return    

def boom_shear_stress_h(dtheta,dx,i,f):#works
    theta = np.linspace(0, 2*math.pi, dtheta)
    shear_flow = np.zeros(dtheta)
    a = np.zeros(dtheta)
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j]+math.pi/2)
    shear_flow = shear_boom(dx,f,i) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress = shear_flow/b_t
    plt.figure(figsize=(19,5))
    plt.suptitle('Shear flow and stress in the mounting boom by horizontal tail')
    plt.subplot(121)
    plt.plot(theta, shear_flow)
    plt.ylabel('Shear flow, N/m')
    plt.xlabel('Location, radian')
    plt.subplot(122)
    plt.plot(theta, shear_stress)
    plt.ylabel('Shear stress, N/m^2')
    plt.xlabel('Location, radian')
    plt.show()
    return

def boom_shear_stress_v(dtheta,dx,i,g):#works
    theta = np.linspace(0, 2*math.pi, dtheta)
    shear_flow = np.zeros(dtheta)
    a = np.zeros(dtheta)
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j])
    shear_flow = force_vert(1000,g) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress = shear_flow/b_t
    plt.figure(figsize=(19,5))
    plt.suptitle('Shear flow and stress in the mounting boom by vertical tail')
    plt.subplot(121)
    plt.plot(theta, shear_flow)
    plt.ylabel('Shear flow, N/m')
    plt.xlabel('Location, radian')
    plt.subplot(122)
    plt.plot(theta, shear_stress)
    plt.ylabel('Shear stress, N/m^2')
    plt.xlabel('Location, radian')
    plt.show()
    return

def torque_vert(dz,g):#works
    z_v = np.linspace(0,v_b,dz)
    v_c = np.linspace(v_r,v_t,dz)
    mid_w = v_w / dz
    w_distr = mid_w * v_c / v_c[dz/2]
    force = g*p.g*w_distr
    moment = force * z_v
    torque = np.sum(moment)
    return torque

def torque_stress(dx, dz, g, i):#works
    t_stress = torque_vert(dz,g) / (2 * b_t * math.pi * r_boom(dx,i)**2)
    return t_stress

def total_shear_stress_boom(dx,dz,dtheta,fh,fv,i):#works
    theta = np.linspace(0, 2*math.pi, dtheta)
    shear_flow = np.zeros(dtheta)
    a = np.zeros(dtheta)
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j])
    shear_flow = shear_boom(dx,fh,i,fv) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress_1 = shear_flow/b_t
    plt.figure(figsize=(19,5))
    plt.suptitle('Total shear stress on the mounting boom')
    plt.subplot(131)
    plt.plot(theta, shear_stress_1)
    plt.ylabel('Shear stress from vertical tail, N/m^2')
    plt.xlabel('Location (from top), radian')
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j]+math.pi/2)
    shear_flow = shear_boom(dx,fv,i,fh) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress_2 = shear_flow/b_t
    plt.subplot(132)
    plt.plot(theta, shear_stress_2)
    plt.ylabel('Shear stress from horizontal tail, N/m^2')
    plt.xlabel('Location (from top), radian')
    shear_stress_tot = shear_stress_1 + shear_stress_2 + torque_stress(dx,dz,fv,i)
    plt.subplot(133)
    plt.plot(theta, shear_stress_tot)
    plt.ylabel('Total shear stress, N/m^2')
    plt.xlabel('Location (from top), radian')
    plt.show()
    return

def total_bending_stress_boom(Ixx,Iyy,dx,dr,i,fh,fv,g):
    x_coord = np.zeros(dr)
    y_coord = np.zeros(dr)
    sin = np.zeros(dr)
    cos = np.zeros(dr)
    theta = np.linspace(0, 2*math.pi, dr)
    for j in range (0,dr):
        sin[j] = math.sin(theta[j]+math.pi/2)
        cos[j] = math.cos(theta[j]+math.pi/2)
    x_coord = r_boom(dx,i) * cos
    y_coord = r_boom(dx,i) * sin
    x = x_boom(dx,i)
    shear = force_vert(1000,fv)
    moment_x = shear * x
    moment_y = moment_boom(dx,fh,g,i)
    stress = np.zeros(dr)
    stress = moment_x/Ixx * y_coord + moment_y/Iyy * x_coord
    plt.figure(figsize=(5,5))
    plt.plot(theta, stress)
    plt.ylabel('Bending stress, N/m^2')
    plt.xlabel('Location (from top), radian')  
    return

def wingbox_MOI(dy,start,end,t,rho):
    
    airfoil_coordinates = np.genfromtxt('foil1_modified.dat',skip_header=1)
    
    Ixx = []
    Iyy = []
    Ixy = []
    
    xcoordinates = np.zeros(len(airfoil_coordinates)) 
    ycoordinates = np.zeros(len(airfoil_coordinates)) 
        
    for i in range(len(airfoil_coordinates)):
        xcoordinates[i] = airfoil_coordinates[i][0]
        ycoordinates[i] = airfoil_coordinates[i][1] 
    
    
    for i in range(0,dy):
        xcoordinates1 = xcoordinates * c_tail(dy,i)
        ycoordinates1 = ycoordinates * c_tail(dy,i)
        
        fit1 = np.polyfit(xcoordinates1[0:134],ycoordinates1[0:134],5)
        f1 = np.poly1d(fit1)
        f11 = f1.deriv(1)
        
        fit2 = np.polyfit(xcoordinates1[133:256],ycoordinates1[133:256],5)
        f2 = np.poly1d(fit2)
        f22 = f2.deriv(1)
    
        front_spar = start * c_tail(dy,i)
        back_spar = end * c_tail(dy,i)
        
        sections = dy            
        dx = 1./sections             
        area_section = dx * t
        x = np.arange(front_spar,back_spar+dx,dx)
    
        arc_length_US = sum(np.sqrt(1+f11(x)**2)*dx) #upper skin
        arc_length_LS = sum(np.sqrt(1+f22(x)**2)*dx) #lower skin
        
        front_spar_length = np.abs(f1(front_spar)) + np.abs(f2(front_spar))
        back_spar_length = np.abs(f1(back_spar)) + np.abs(f2(back_spar))
        
        front_spar_area = np.abs(f1(front_spar)) + np.abs(f2(front_spar))*t
        back_spar_area = np.abs(f1(back_spar)) + np.abs(f2(back_spar))*t
        
        total_area = (arc_length_US + arc_length_LS)*t + front_spar_area + back_spar_area
        
        y_NA = (area_section * (sum(f1(x)) + sum(f2(x))) + front_spar_area*(f1(front_spar) + f2(front_spar))/2 + back_spar_area*(f1(back_spar) + f2(back_spar))/2)/ total_area 
        x_NA = (area_section * sum(x) * 2 + front_spar_area * front_spar + back_spar_area*back_spar) / total_area
                     
        Ixx_section = 2*dx*t**3/12 + area_section * sum(f1(x)**2 + f2(x)**2) + t/12 * (front_spar_length**3 + front_spar_area*(f1(front_spar) + f2(front_spar - y_NA))**2 + back_spar_length**3 + back_spar_area*(f1(back_spar) + f2(back_spar - y_NA))**2)
        Iyy_section = (front_spar_area * (front_spar-x_NA)**2) + (back_spar_area * (back_spar-x_NA)**2) + area_section*(sum(x**2)-2*x_NA*sum(x)+len(x)*x_NA**2)
        Ixy_section = area_section*(sum((x-x_NA)*(f1(x)-y_NA)) + sum((x-x_NA)*(f2(x)-y_NA))) + front_spar_area * (front_spar - x_NA) * ((f1(front_spar) + f2(front_spar))/2 - y_NA) + back_spar_area * (back_spar - x_NA) * ((f1(back_spar) + f2(back_spar))/2 - y_NA)
    
        Ixx.append(Ixx_section)
        Iyy.append(Iyy_section)
        Ixy.append(Ixy_section)
        
    weight = np.sum(total_area*rho*p.b_ht/dy)
        
    return Ixx, Iyy, Ixy, weight, f1, f2, y_NA, x_NA, x 
#                master_function(500,500,100,100,100,100,0.15,0.575,0.002,0,0)

def master_function(fh,fv,dx,dy,dz,dtheta,start,end,t,m1,m2):
    #inputs: wingbox material, boom material numbers
    Ixx = wingbox_MOI(dy,start,end,t,mat.rho[m1])[0]
    Ixy = wingbox_MOI(dy,start,end,t,mat.rho[m1])[2]
    Iyy = wingbox_MOI(dy,start,end,t,mat.rho[m1])[1]
    tail_plots(dy,fh)
    tail_shear_stress(dy,dx,dz,dy-1,fh,Ixx,Ixy,Iyy,t)
    boom_plots(dx,fh,1)
    boom_plots_vertical(dx,fv)
    total_shear_stress_boom(dx,dz,dtheta,fh,fv,dz-1)
    total_bending_stress_boom(Ixx,Iyy,dx,dtheta,dx-1,fh,fv,1)
    print ('Weight of the wingbox:')
    print wingbox_MOI(dy,start,end,t,mat.rho[m1])[3]
    print ('Weight of the boom:')
    print W_boom(mat.rho[m2])
    return