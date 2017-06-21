import numpy as np

import sys
sys.path.insert(0,'C:\Users\Claudia\Models')
import Init_Parameters as p
import Wing_model as w
import Fuselage_composite_final as fc

## FUSELAGE PARAMETERS
L_f     = p.Lf   # m 
                 # This depends on the aerodynamic mesh and the location of the
                 # forces and moments
c_r    = w.C1
n       =  2*6


## LOADS ON FUSELAGE
dp      = 101325-50600.9 # Pressure difference at 18000 ft
p_fus   = 0.1985
x_cg    = 0.42*p.MAC+ p.X_le_mac
x_ac    = 0.25*p.MAC+ p.X_le_mac
q_fus   = n*((p.MTOW*9.81)*p_fus)/L_f

L_w_y   = w.wing_shear(w.CL,p.rho_0,p.V_cruise)[-1][0]
L_h_y   = q_fus*L_f-L_w_y

q_w     = 3*L_w_y/c_r
L_h_y   = q_fus*L_f-L_w_y
x_ac_back  = L_f-(x_ac + c_r*2./3)
x_ac_front = L_f-(x_ac - c_r*1./3)

# MESH OF THE FUSELAGE @ EVERY CM ALONG FUSELAGE LENGTH 
dl      = 0.01   # cm
# For ease of calculation all sections are assumed to taper uniformly

## SECTION 1: ENGINE MOUNT + FIREWALL:
r_i       = 0.4    #m
L_i       = 1.35   #m

x_i       = np.linspace(0,L_i,L_i/dl)
r_s_i     = r_i*np.ones(np.shape(x_i)[0])

## SECTION 2: COCKPIT
r_ii_front = r_i      #m
r_ii_max   = 0.7      #m
r_ii_back  = 0.4      #m

L_ii = 2.7            #m         
L_ii_front    = 1.3   #m
L_ii_const    = 0.9   #m
L_ii_back     = L_ii-L_ii_front-L_ii_const #m 


x_front   = np.linspace(L_i,L_ii_front+L_i,L_ii_front/dl)
x_const   = np.linspace(L_ii_front+L_i,L_ii_front+L_i+L_ii_const,L_ii_const/dl)
x_back    = np.linspace(L_ii_front+L_i+L_ii_const,L_ii_front+L_i+ \
                        L_ii_const+L_ii_back,L_ii_back/dl)

r_s_ii_front = np.linspace(r_ii_front,r_ii_max,np.shape(x_front)[0])
r_s_ii_const = r_ii_max*np.ones(np.shape(x_const)[0])
r_s_ii_back  = np.linspace(r_ii_max,r_ii_back,np.shape(x_back)[0])

## SECTION 3: BACK PART
r_iii_front       = r_ii_back #m
r_iii_back        = 0.1

L_iii             = L_f -L_i-L_ii

x_iii             = np.linspace(L_ii+L_i,L_f,L_iii/dl)
r_s_iii           = np.linspace(r_iii_front,r_iii_back,np.shape(x_iii)[0])

## COMBINED  SECTIONS 1+2+3
x = np.concatenate((x_i,x_front))
x = np.append(x,x_const)
x = np.append(x,x_back)
x = np.append(x,x_iii)


r= np.concatenate((r_s_i,r_s_ii_front))
r= np.append(r,r_s_ii_const)
r= np.append(r,r_s_ii_back)
r= np.append(r,r_s_iii)

##S_y = np.zeros(np.shape(x))
##M_x = np.zeros(np.shape(x))
##for space in range(np.shape(x)[0]):
##    if x[space] < (L_f-x_ac):
##        S_y[space] = q_fus*x[space]-L_h_y
##        M_x[space] = q_fus*x[space]**2/2-L_h_y*x[space]
##    if x[space] > (L_f-x_ac):
##        S_y[space] = q_fus*x[space]-L_w_y-L_h_y
##        M_x[space] = q_fus*x[space]**2/2-L_w_y*(x[space]-x_ac)-L_h_y*x[space]
##S_x = np.zeros(np.shape(x)[0])
##M_y = np.zeros(np.shape(x)[0])
##T= np.zeros(np.shape(x)[0])

S_y = np.zeros(np.shape(x))
M_x = np.zeros(np.shape(x))
for space in range(np.shape(x)[0]):
    if x[space] < x_ac_back:
        S_y[space] = q_fus*x[space]-L_h_y
        M_x[space] = q_fus*x[space]**2/2-L_h_y*x[space]
    if x[space] > x_ac_back and x[space] < x_ac_front:
        S_y[space] = q_fus*x[space]-q_w*c_r/2*(1/3-(x[space]-x_ac_back)**2/ \
                                               c_r**2)-L_h_y
        M_x[space] = q_fus*x[space]**2/2-q_w*c_r*(x[space]-x_ac_back)/6* \
                     (1-(x[space]-x_ac_back)**2/c_r**2)-L_h_y*x[space]
    if x[space] > x_ac_front:
        S_y[space] = q_fus*x[space]-L_w_y-L_h_y
        M_x[space] = q_fus*x[space]**2/2-L_w_y*(x[space]-x_ac_back)- \
                     L_h_y*x[space]
M_y = np.zeros(np.shape(x)[0])
S_x = np.zeros(np.shape(x)[0])
T= np.zeros(np.shape(x)[0])


## Material Properties Core Material
E_c     = 400*10**6
t_c_min = 3*10**(-3)
G_xz    = 85*10**6
G_yz    = 50*10**6
rho_c   = 96.

## Material Properties Facing Material
E_x     = 60.1*10**9
E_y     = 60.1*10**9
v_xy    = 0.307
rho_f   = 1580.
sigma_t = 365.*10**6
sigma_c = 657.*10**6
e_t     = sigma_t/E_x
e_c     = sigma_c/E_y
t_f_min = 0.3
tau     = (-sigma_t+sigma_c)/2
sigma_hoop    = 100*10**6
d =fc.thickness(r,M_x,M_y,S_x,S_y,T,sigma_t,sigma_c,\
    tau,t_c_min,dp,sigma_hoop,rho_f,rho_c,E_x,E_c)
m = np.sum(fc.mass(r,M_x,M_y,S_x,S_y,T,sigma_t,sigma_c,tau,\
                   t_c_min,dp,sigma_hoop,rho_f,rho_c,E_x,E_c))
print m   
