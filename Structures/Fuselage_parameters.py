import numpy as np
#import Init_parameters as i
## FUSELAGE PARAMETERS
t_min   = 0.5    #mm
L_f     = 9.04   # m 
                 # This depends on the aerodynamic mesh and the location of the
                 # forces and moments
dl      = 0.01   # mm
dp      = 101325-50600.9 # Pressure difference at 18000 ft
#W       = i.MTOW*9.81    # Fully Loaded Aircraft
# For ease of calculation all sections are assumed to taper uniformly 
## SECTION 1: ENGINE MOUNT + FIREWALL:
r_i       = 0.4    #m
L_i       = 1.35   #m
x_i       = np.linspace(0,L_i,L_i/dl)
r_s_i       = r_i*np.ones(np.shape(x_i)[0])
## SECTION 2: COCKPIT
r_ii_front = r_i      #m
r_ii_max   = 0.7      #m
r_ii_back  = 0.4      #m
L_ii = 3.0            #m         
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
r_iii_front       = r_ii_back    #m
r_iii_back        = 0.1   
L_iii             = L_f -L_i-L_ii   #m
x_iii             = np.linspace(L_ii+L_i,L_f,L_iii/dl)
r_s_iii           = np.linspace(r_iii_front,r_iii_back,np.shape(x_iii)[0])

## 
x = np.concatenate((x_i,x_front))
x = np.append(x,x_const)
x = np.append(x,x_back)
x = np.append(x,x_iii)
##
r= np.concatenate((r_s_i,r_s_ii_front))
r= np.append(r,r_s_ii_const)
r= np.append(r,r_s_ii_back)
r= np.append(r,r_s_iii)

## Material Properties Core Material
E_c     = 400*10**6
t_c_min = 3*10**(-3)
G_xz    = 85*10**6
G_yz    = 50*10**6
rho_c   = 96

## Material Properties Facing Material
E_x     = 60.1*10**9
E_y     = 60.1*10**9
v_xy    = 0.307
rho_f   = 1580
sigma_t = 356*10**6
sigma_c = 657*10**6
e_t     = sigma_t/E_x
e_c     = sigma_c/E_y
t_f_min = 0.3
tau     = 11.5*10**6
hoop    = 100*10**6

