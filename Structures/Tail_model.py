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

"""


"""
Stuff

"""
def y_tail(dy,i):
    y_tail = np.linspace(0,p.b_ht/2,dy)
    return y_tail[i]

def y_tail_inv(dy,i):
    y_tail_inv = np.linspace(p.b_ht/2,0,dy)
    return y_tail_inv[i]

def c_tail(dy,i):
    c_tail = np.linspace(p.ct_ht,p.cr_ht,dy+1)
    return c_tail[i]

def w_tail(dy,i):
    mid_w = p.W_ht / dy
    w_distr = mid_w * c_tail(dy,i) / c_tail(dy,dy/2)
    return w_distr

def slope_c(dy):
    slope = (p.ct_ht+p.cr_ht)/dy
    return slope

def tail_I_xx(dy,i):
    I_xx = np.zeros(dy)
    I_xx = (p.cr_ht - y_tail_inv(dy,i) * slope_c(dy))**4 * (box_h)**3 * box_b/12 - (p.cr_ht - y_tail_inv(dy,i) * slope_c(dy))**4 * (box_h - box_t)**3 * (box_b - 0.0005)/12
    return I_xx
    #9G ULTIMATE STATIC LOAD HORIZONTAL TAIL

#Distribution at every point
def tail_force(dy,g,i):
    force = g*p.g*w_tail(dy,i)
    return force

def tail_distr(dy,g,i):
    distr = g*p.g*w_tail(dy,i)/(p.b_ht/2*dy)
    return distr 

def tail_shear(dy,g,j):
    shear = np.zeros(j+1)
    shear[0] = tail_force(dy,g,0)
    for i in range (1,j+1):
        shear[i] = shear[i-1] + tail_force(dy,g,i)
    return shear[j]

def tail_moment(dy,g,j):
    moment = np.zeros(j+1)
    moment[0] = tail_force(dy,g,0) * y_tail(dy,0)
    for i in range (1,j+1):
        moment[i] = moment[i-1] + tail_force(dy,g,i) * y_tail(dy,i)
    return moment[j]
        
def tail_plots(dy,g):
    distr = np.zeros(dy)
    y = np.zeros(dy)
    shear = np.zeros(dy)
    moment = np.zeros(dy)
    for i in range (0,dy):
        distr[i] = tail_distr(dy,g,i)
    for i in range (0,dy):
        y[i] = y_tail(dy,i)
    for i in range (0,dy):
        shear[i] = tail_shear(dy,g,i)
    for i in range (0,dy):
        moment[i] = tail_moment(dy,g,i)
    plt.figure(figsize=(19,5))
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

def tail_shear_stress(dy,dx,dz,i,g):
    x_12 = np.linspace(0,((p.cr_ht - y_tail_inv(dy,i) * slope_c(dy)) * box_b/2),dx)
    z_23 = np.linspace(((p.cr_ht - y_tail_inv(dy,i) * slope_c(dy)) * box_h)/2, -1*(p.cr_ht-y_tail_inv(dy,i)*slope_c(dy))*0.1*p.cr_ht/2, dz)
    x_34 = np.linspace((p.cr_ht - y_tail_inv(dy,i) * slope_c(dy)) * box_b/2,0,dx)
    q_12 = np.zeros(dx)
    q_23 = np.zeros(dy)
    q_34 = np.zeros(dx)
    for j in range (0,dx):
        q_12[j] = -tail_shear(dy,g,i) * box_t * z_23[0] * x_12[j] / (2*tail_I_xx(dy,i))
    for j in range (0,dz):
        q_23[j] = q_12[dx-1] + box_t * tail_shear(dy,g,i) * x_12[dx-1] * ((z_23[0])**2 - z_23[j]**2)/(4*tail_I_xx(dy,i))
    for j in range (0,dx):
        q_34[j] = q_23[dz-1] - tail_shear(dy,g,i) * box_t * z_23[dz-1] * x_12[j] / (2*tail_I_xx(dy,i))
    tau_12 = np.zeros(dx)
    tau_23 = np.zeros(dx)
    tau_34 = np.zeros(dx)
    for j in range (0,dx):
        tau_12[j] = q_12[j] / box_t
    for j in range (0,dx):
        tau_23[j] = q_23[j] / box_t
    for j in range (0,dx):
        tau_34[j] = q_34[j] / box_t
    plt.figure(figsize=(19,5))
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

def y_boom(dx,i):
    y = np.linspace(0,b_l,dx)
    return y[i]

def r_boom(dx,i):
    slope = (b_rr-b_rt)/dx
    r = b_rt + slope*i
    return r

def W_boom(rho):
    w = rho * (math.pi * (b_rr+b_rt)**2/4 * b_l - math.pi * (b_rr+b_rt-2*b_t)**2/4*b_l)
    return w

def w_boom(dx,i):
    mid_w = W_boom(mat.rho[0]) / dx
    w_distr = mid_w * r_boom(dx,i) / r_boom(dx,dx/2)
    return w_distr

def force_boom(dx,g,i):
    force = g*p.g*w_boom(dx,i)
    return force

def distr_boom(dx,g,i):
    distr = (g*p.g*w_boom(dx,i) + tail_shear(1000,g,999))/(b_l*dx) 
    return distr

def shear_boom(dx,g,j):
    distr = np.zeros(j+1)
    for i in range (1,j+1):
        distr[i] = distr_boom(dx,g,i) 
    shear_j = force_boom(dx,g,0) * (j+1) + tail_shear(1000,g,999) + sum(distr)
    return shear_j

def moment_boom(dx,g,j):
    moment = np.zeros(j+1)
    moment[0] = (force_boom(dx,g,0) + tail_shear(1000,g,999)) * y_tail(dx,0)
    for i in range (1,j+1):
        moment[i] = moment[i-1] + tail_force(dx,g,i) * y_tail(dx,i)
    return moment[j]
        
def boom_plots(dx,g):
    distr = np.zeros(dx)
    y = np.zeros(dx)
    shear = np.zeros(dx)
    moment = np.zeros(dx)
    for i in range (0,dx):
        distr[i] = distr_boom(dx,g,i)
    for i in range (0,dx):
        y[i] = y_boom(dx,i)
    shear[0] = shear_boom(dx,g,0)
    for i in range (1,dx):
        shear[i] = shear[i-1] + distr_boom(dx,g,i)
    for i in range (0,dx):
        moment[i] = moment_boom(dx,g,i)
    plt.figure(figsize=(19,5))
    plt.subplot(131)
    plt.plot(y, distr)
    plt.ylabel('Load, N')
    plt.xlabel('Location, m')
    plt.subplot(132)
    plt.plot(y, shear)
    plt.ylabel('Shear, N')
    plt.xlabel('Location, m')
    plt.subplot(133)
    plt.plot(y, moment)
    plt.ylabel('Moment, Nm')
    plt.xlabel('Location, m')
    plt.show()
    return    
    
def boom_shear_stress_h(dtheta,dx,i,g):
    theta = np.linspace(0, 2*math.pi, dtheta)
    shear_flow = np.zeros(dtheta)
    a = np.zeros(dtheta)
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j]+math.pi/2)
    shear_flow = shear_boom(dx,g,i) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress = shear_flow/b_t
    plt.figure(figsize=(19,5))
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

def force_vert(dz,g):
    v_c = np.linspace(v_r,v_t,dz)
    mid_w = v_w / dz
    w_distr = mid_w * v_c / v_c[dz/2]
    force = g*p.g*w_distr
    return sum(force)

def boom_shear_stress_v(dtheta,dx,i,g):
    theta = np.linspace(0, 2*math.pi, dtheta)
    shear_flow = np.zeros(dtheta)
    a = np.zeros(dtheta)
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j])
    shear_flow = shear_boom(dx,g,i) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress = shear_flow/b_t
    plt.figure(figsize=(19,5))
    plt.subplot(121)
    plt.plot(theta, shear_flow)
    plt.ylabel('Shear flow, N/m')
    plt.xlabel('Location, radian')
    plt.subplot(122)
    plt.plot(theta, shear_stress)
    plt.ylabel('Shear stress, N/m^2')
    plt.xlabel('Location, radian')
    plt.show()

def torque_vert(dz,g):
    z_v = np.linspace(0,v_b,dz)
    v_c = np.linspace(v_r,v_t,dz)
    mid_w = v_w / dz
    w_distr = mid_w * v_c / v_c[dz/2]
    force = g*p.g*w_distr
    moment = force * z_v
    torque = sum(moment)
    return torque

def torque_stress(dx, dz, g, i):
    t_stress = torque_vert(dz,g) / (2 * b_t * math.pi * r_boom(dx,i)**2)
    return t_stress

def total_shear_stress_boom(dx,dz,dtheta,g1,g2,i):
    theta = np.linspace(0, 2*math.pi, dtheta)
    shear_flow = np.zeros(dtheta)
    a = np.zeros(dtheta)
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j])
    shear_flow = shear_boom(dx,g1,i) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress_1 = shear_flow/b_t
    plt.figure(figsize=(19,5))
    plt.subplot(131)
    plt.plot(theta, shear_stress_1)
    plt.ylabel('Shear stress from vertical tail, N/m^2')
    plt.xlabel('Location (from top), radian')
    for j in range (0,dtheta):
        a[j] = math.cos(theta[j]+math.pi/2)
    shear_flow = shear_boom(dx,g2,i) * (a-1) / (math.pi * r_boom(dx,i))
    shear_stress_2 = shear_flow/b_t
    plt.subplot(132)
    plt.plot(theta, shear_stress_2)
    plt.ylabel('Shear stress from horizontal tail, N/m^2')
    plt.xlabel('Location (from top), radian')
    print torque_stress(dz,g2)
    shear_stress_tot = shear_stress_1 + shear_stress_2 + torque_stress(dx,dz,g2,i)
    plt.subplot(133)
    plt.plot(theta, shear_stress_tot)
    plt.ylabel('Total shear stress, N/m^2')
    plt.xlabel('Location (from top), radian')
    plt.show()

def big_scary_funtion():
        
    






#Deflection
#a = 0.5 * forces_9g[999]
#b = 0.25 * (distr_9g[999] - distr_9g[998])/p.b_ht/2000
#c = p.cr_ht
#d = slope_c
#e = mat.E[0] * (0.1*p.cr_ht)**3 * 0.25*p.cr_ht/12 - (0.1*p.cr_ht - 0.0005)**3 * (0.25*p.cr_ht/12 - 0.0005)
#v_9g = (-(b*d*y_tail+a*d-4*b*c)*np.log(abs(d*y_tail-c)))/(e*d*d*d*d*d) + (c*(3*(2*a*d-3*b*c)*d*y_tail-(5*a*d-8*b*c)*c))/(6*d*d*d*d*d*e*(d*y_tail-c)*(d*y_tail-c)) + (b*y_tail)/(d*d*d*d*e) + y_tail*(b*np.log(abs(-c)))/(e*d*d*d*d) - y_tail*(2*a*d-11*b*c)*c*c/(6*d*d*d*d*c*c*c*e) + (a*d-4*b*c)*np.log(abs(-c))/(d*d*d*d*d*e) + ((5*a*d-8*b*c)*c)/(6*d*d*d*d*d*c*c*e)       
#plt.subplot(144)
#plt.plot(y_tail, v_9g)
#plt.show()       
#Deflection
#deflec_45g = forces_45g * (p.b_ht/4000)**2 * (3 * p.b_ht/2000-p.b_ht/4000) /(6 * mat.E[0] * I_xx)
#for i in range (2,dy+1):
#    deflec_45g[dy-i] = deflec_45g[dy-i] + deflec_45g[dy-i+1]
#plt.subplot(144)
#plt.plot(y_tail, deflec_45g)
#plt.show()
