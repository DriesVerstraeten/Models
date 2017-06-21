# -*- coding: utf-8 -*-
"""
Created on Wed Jun 07 14:17:58 2017

@author: Maxime
"""

import numpy as np
import matplotlib.pyplot as plt

## Create folder to store plots in       
#import os
#directory = os.getcwd() + '\Plots'
#if not os.path.exists(directory):
#    os.makedirs(directory)

#### INPUTS ##########
#    OEW = 997. # [kg]
#    W_wing = wing mass [kg]
#    c = 1.568 # MAC
#    x_ac = 0.3
#    downwash = 0.427
#    Cmac = -0.121
#    CL_Ah_land = 2.3
#    CL_alpha_Ah_land = 5.11
#    CL_h_land = -0.55
#    CL_alpha_h_land = 3.852
#    VhV2 = 0.85**2 # (Vh/V)^2
#    xcg_OEW
#    xcg_wing
#    lh = tail arm [m]

##### OUTPUTS StabTailSizing() ####
#    x_lemac_array = array with possible x_lemac locations [m]
#    xcg_for_load = xcg_locations during loading wrt to MAC for different x_LEMAC
#    xcg_aft_load
#    xcg_pil_load
#    mass_tot = total mass during loading [kg]
#    mass_tot_pil
#    xcg_for = most forward xcg locations wrt MAC (incl 2% MAC margin)
#    xcg_aft
#    xcg_control = xcg/MAC control curve for different Sh/S
#    xcg_stab = xcg/MAC stability curve for different Sh/S
#    ShS_array = array with possible Sh/S ratios
#    ShS_opt = optimal Sh/S
#    x_lemac_opt = optimal x_LEMAC [m]
#    xcg_for_opt = most forward xcg/MAC for optimal x_LEMAC configuration
#    xcg_aft_opt 
#    where_opt = index of optimal configuration
#################



def StabTailSizing(OEW,W_wing,c,x_ac,downwash,Cmac,CL_Ah_land,CL_alpha_Ah_land,CL_h_land,CL_alpha_h_land,VhV2,xcg_OEW,xcg_wing,lh, x_lemac_acc, ShS_acc): # Calc optimal Sh/S and x_lemac. Input accuracy of analysis.
    ##### Calc xcg_OEW for different x_lemacs ###
#    xcg_OEW = 2.9 # [m] # Starting value. Must be updated with xcg_OEW_opt
#    xcg_wing = 2.9 # [m] # Starting value 
    xcg_wingMAC = 0.3 # cg of wing wrt MAC
#    OEW = 997. # [kg]
#    W_wing = 209. # [kg]
    xcg_fuelMAC = 0.35 # xcg of fuel wrt MAC
#    c = 1.568 # MAC
    
    moment_Aw = xcg_OEW * OEW - xcg_wing * W_wing    
    x_lemac_array = np.arange(1.5,3.0,x_lemac_acc) # Set range and accuracy of analysis !!!
    xcg_fuel_array = x_lemac_array + xcg_fuelMAC * c
    xcg_OEW_array = (moment_Aw + (x_lemac_array + xcg_wingMAC * c) * W_wing) / OEW    
    ###########################
    ##### Calc most forward and most aft xcg for different loading situations ###
    xcg_cargo = 3.75 # [m]
    xcg_pilot = 2.57
    xcg_pax = 3.51

    cargo = 100. # [kg]
    pilot = 86.
    pax = 86.
    fuel = 219.
    
    xcg = [xcg_cargo, xcg_pilot, xcg_pax, xcg_fuel_array]
    mass = [cargo, pilot, pax, fuel]
    
    order_for = [0,1,1,2,2,3] # Full loading. Pilots first    
    xcg_for_load = [(xcg_OEW_array - x_lemac_array) / c]
    mass_tot = [OEW]
    moment = xcg_OEW_array * OEW
    for i in order_for:
        moment += xcg[i] * mass[i]
        mass_tot.append(mass_tot[-1] + mass[i])
        xcg_for_load.append(((moment / mass_tot[-1])-x_lemac_array) / c)
        
    order_aft = [0,2,2,1,1,3] # Full loading. Passengers first   
    xcg_aft_load = [(xcg_OEW_array - x_lemac_array) / c]
    mass_tot = [OEW]
    moment = xcg_OEW_array * OEW
    for i in order_aft:
        moment += xcg[i] * mass[i]
        mass_tot.append(mass_tot[-1] + mass[i])
        xcg_aft_load.append(((moment / mass_tot[-1])-x_lemac_array) / c)
    
    order_pilotsonly = [1,1,3] # Light loading. Only pilots. No passengers. No cargo 
    xcg_pil_load = [(xcg_OEW_array - x_lemac_array) / c]
    mass_tot_pil = [OEW]
    moment = xcg_OEW_array * OEW
    for i in order_pilotsonly:
        moment += xcg[i] * mass[i]
        mass_tot_pil.append(mass_tot_pil[-1] + mass[i])
        xcg_pil_load.append(((moment / mass_tot_pil[-1])-x_lemac_array) / c)
    # Find most forward and most aft position for every x_lemac
    xcg_for = np.amin(xcg_for_load + xcg_aft_load + xcg_pil_load, axis=0) - 0.02
    xcg_aft = np.amax(xcg_for_load + xcg_aft_load + xcg_pil_load, axis=0) + 0.02
    #############################
    ##### Calc control and stability curves ###
#    x_ac = 0.3
#    downwash = 0.427
#    Cmac = -0.121
#    CL_Ah_land = 2.3
#    CL_alpha_Ah_land = 5.11
#    CL_h_land = -0.55
#    CL_alpha_h_land = 3.852
#    lh = 5.425 # Tail arm
#    VhV2 = 0.85**2 # (Vh/V)^2
    SM = 0.05 # Stability Margin
 
    ShS_array = np.arange(0,0.4,ShS_acc) # Sh/S. Set range and accuracy of analysis !!!
    
    xcg_control = x_ac - Cmac/CL_Ah_land + CL_h_land/CL_Ah_land * ShS_array * lh / c * VhV2
    xcg_stab = x_ac + CL_alpha_h_land / CL_alpha_Ah_land * (1-downwash) * ShS_array * lh / c * VhV2 - SM
    #############################
    ##### Determine minimum Sh/S for optimal x_lemac ###    
    ShS_min_array = []
    ShS_diff_array = []
    for i in range(len(xcg_for)):
        where_for = np.where(xcg_control < xcg_for[i])[0] # Loc of intersect of control curve with most forward xcg
        where_aft = np.where(xcg_stab > xcg_aft[i])[0] # Loc of intersect of stability curve with most aft xcg
        if len(where_for) == 0 or len(where_aft) == 0:
            ShS_min_array.append(9)
            ShS_diff_array.append(9)
            continue
        con, stab = ShS_array[where_for[0]], ShS_array[where_aft[0]] # Sh/S value of intersections
        ShS_min_array.append(max(con, stab))
        ShS_diff_array.append(abs(con - stab))
    where_opt = np.argmin(ShS_diff_array)
    ShS_opt = ShS_min_array[where_opt]
    x_lemac_opt = x_lemac_array[where_opt]
    xcg_OEW_opt = xcg_OEW_array[where_opt]
    xcg_for_opt, xcg_aft_opt = xcg_for[where_opt], xcg_aft[where_opt]
    return x_lemac_array, xcg_for_load, xcg_aft_load, xcg_pil_load, mass_tot, mass_tot_pil, xcg_for, xcg_aft, xcg_control, xcg_stab, ShS_array, ShS_opt, x_lemac_opt, xcg_OEW_opt, xcg_for_opt, xcg_aft_opt, where_opt 



def StabTailSizingPlots(OEW,W_wing,c,x_ac,downwash,Cmac,CL_Ah_land,CL_alpha_Ah_land,CL_h_land,CL_alpha_h_land,VhV2,xcg_OEW,xcg_wing,lh, x_lemac_acc, ShS_acc):
    # Create folder to store plots in       
    import os
    directory = os.getcwd() + '\Plots'
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    x_lemac_array, xcg_for_load, xcg_aft_load, xcg_pil_load, mass_tot, mass_tot_pil, xcg_for, xcg_aft, xcg_control, xcg_stab, ShS_array, ShS_opt, x_lemac_opt, xcg_OEW_opt, xcg_for_opt, xcg_aft_opt, where_opt = StabTailSizing(OEW,W_wing,c,x_ac,downwash,Cmac,CL_Ah_land,CL_alpha_Ah_land,CL_h_land,CL_alpha_h_land,VhV2,xcg_OEW,xcg_wing,lh, x_lemac_acc, ShS_acc)
    
    ##### Plot Loading diagram with all x_lemac options ###
#    plt.plot(xcg_for_load,mass_tot, xcg_aft_load,mass_tot, xcg_pil_load,mass_tot_pil)
#    plt.xlabel('x_cg/MAC')
#    plt.ylabel('Total Mass [kg]')
#    plt.savefig(directory + '\Loading.png', format='png', bbox_inches='tight', dpi=400)
#    plt.show()
    ########################
    ##### Plot Loading diagram for optimal x_lemac ###
    fig = plt.plot(np.array(xcg_for_load)[:,where_opt],mass_tot, np.array(xcg_aft_load)[:,where_opt],mass_tot, np.array(xcg_pil_load)[:,where_opt],mass_tot_pil)
    plt.xlabel('x_cg/MAC')
    plt.ylabel('Total Mass [kg]')
    plt.legend(iter(fig), ('Full Loading (Pilots first)', 'Full Loading (Pax first))', 'Light Loading (2 pilots only)'), bbox_to_anchor=(1.7,1.0)) #loc=2)
    plt.savefig(directory + '\LoadingOpt.png', format='png', bbox_inches='tight', dpi=400)
    plt.show()
    #######################
#    plt.plot(xcg_control,ShS_array, xcg_stab,ShS_array)
#    plt.xlabel('x_cg/MAC')
#    plt.ylabel('Sh/S')
#    plt.savefig(directory + '\StabControl.png', format='png', bbox_inches='tight', dpi=400)
#    plt.show()
    
    ##### Plot stability and control curves with x_lemac options ###
    fig, ax1 = plt.subplots(figsize=(10,5))
    ax1.plot(xcg_control, ShS_array, color='green', label='Controllability')
    ax1.plot(xcg_stab, ShS_array, color='green', linestyle='dashed', label='Stability')
    ax1.plot([0,1],[ShS_opt,ShS_opt],color='black',linestyle='dotted')
    opt_ratio = ShS_opt / (ShS_array[-1] - ShS_array[0])
    ax2lim = (x_lemac_opt - x_lemac_array[0]) / opt_ratio + x_lemac_array[0]
    ax1.set_xlim(0,1)
    ax1.set_xlabel('x_cg/MAC')
    ax1.set_ylabel('Sh/S', color='g')
    ax1.tick_params('y', colors='g')
    plt.legend(loc=2)
    
    ax2 = ax1.twinx()
    ax2.plot(xcg_for,x_lemac_array,color='red',label='Most forward x_cg')
    ax2.plot(xcg_aft,x_lemac_array,color='red',linestyle='dashed', label='Most aft x_cg')
    ax2.set_xlim(0,1)
    ax2.set_ylim([x_lemac_array[0], ax2lim])
    ax2.set_ylabel('x_LEMAC [m]', color='r')
    ax2.tick_params('y', colors='r')
    plt.legend()
    plt.savefig(directory + '\StabControl+x_lemac.png', format='png', bbox_inches='tight', dpi=400)
    plt.show() 
    return



#### RUN FUNCTIONS: ####
#StabTailSizingPlots(997,209,1.568,0.3,0.427,-0.121,2.3,5.11,-0.55,3.852,0.85**2,2.9,2.9,5.425,0.001,0.001)
#
x_lemac_array, xcg_for_load, xcg_aft_load, xcg_pil_load, mass_tot, mass_tot_pil, xcg_for, xcg_aft, xcg_control, xcg_stab, ShS_array, ShS_opt, x_lemac_opt, xcg_OEW_opt, xcg_for_opt, xcg_aft_opt, where_opt = StabTailSizing(997,209,1.568,0.3,0.427,-0.121,2.3,5.11,-0.55,3.852,0.85**2,2.9,2.9,5.425,0.001,0.001)
############