import numpy as np
## FUSELAGE PARAMETERS
L_f     = 9.04   # m 
dl      = 0.01   # Spacing between sections in fuselage
                 # This depends on the aerodynamic mesh and the location of the
                 # forces and moments
## SECTION 1: ENGINE MOUNT + FIREWALL: For this section only the width will
                 # change. 
a_i_front = 0.4  #m 
a_i_end   = 0.6  #m
b       = 0.3    #m 
L_i     = 1.35   #m
a_i = np.linspace(a_i_front,a_i_end,L_i/dl) # 
b_i = b*np.ones(np.shape(a_i))              #

## SECTION 2: COCKPIT

a_ii = a_i_end       #m
b_ellipse_top = b    #m
b_ellipse_bot = b    #m
L_ii = 3.0           #m      
D    = 0.8           #m    
L_ii_front    = 1.3  #m
L_ii_const    = 0.9  #m 
L_ii_back     = L_ii-L_ii_front-L_ii_const #m
d_front   = np.linspace(0,D,L_ii_front/dl)
d_const   = D*np.ones(L_ii_const/dl)
d_back    = np.linspace(D,0,L_ii_back/dl)

## SECTION 3: BACK PART

a_iii_front = a_ii            #m
b_iii_front = b_ellipse_bot   #m
a_iii_end   = 0.1             #m
b_iii_end   = 0.1             #m
L_iii = L_f-L_i-L_ii          #m
a_iii = np.linspace(a_iii_front,a_iii_end,L_iii/dl)
b_iii = np.linspace(b_iii_front,b_iii_end,L_iii/dl)
