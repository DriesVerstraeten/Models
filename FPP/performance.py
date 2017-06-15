# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 17:26:23 2017

@author: Simon

This is the flight performance file where all functions come together

Basically, if one needs power calculations, the propulsion file will give them.

The flight performance is divided into cruise perfornmance and climb, descent 
and field performance.
"""

import propulsion as prop
import Climb_Descent_FieldPerformance as CDFP
import CruiseRange as CR

def Performance(Mfuelmax, PL, V_acc, MTOW, OEW, S, A, e, CD_0, etha_p, etha_p_landing, P_TO, P, T_static, LD_max, hcruise, CL_0, CL_alpha, alpha_TO, CL_max_TO, CL_landing_max, CL_landing_td, D_prop, D_spinner, V_c, V_H):
    
    R_eff, R_spec, R_Vmax, V_Cruise_eff, V_Cruise_spec, V_max, PowerReq, PowerAva, PowerCruise_spec, V_min, V_Loiter, V_array = CR.calcCruise(MTOW,OEW,Mfuelmax,PL,S,V_acc)
    CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si, t_cruiseh, ROC_SL_si, h_absolute_si, R_glide, S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO = CDFP.main_updown(MTOW, S, A, e, CD_0, etha_p, etha_p_landing, P_TO, P, T_static, LD_max, hcruise, CL_0, CL_alpha, alpha_TO, CL_max_TO, CL_landing_max, CL_landing_td, D_prop, D_spinner, V_c, V_H)
    
    return 


"""
def main_updown(MTOW, S, A, e, CD_0, etha_p, etha_p_landing, P_TO, P, T_static, LD_max, hcruise, CL_0, CL_alpha, alpha_TO, CL_max_TO, CL_landing_max, CL_landing_td, D_prop, D_spinner, V_c, V_H):
    CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si = Takeoff(W_TO, S, A, e, CD_0, etha_p, P_TO, CL_0, CL_alpha, alpha_TO, CL_max_TO, D_prop, D_spinner, V_c, V_H)
    t_cruiseh, ROC_SL_si, h_absolute_si = Climb(W_TO, S, A, e, CD_0, etha_p, P, LD_max, hcruise)
    R_glide = Descent(W_TO, S, A, e, CD_0, LD_max, hcruise)
    S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO = Landing(W_TO, S, A, e, CD_0, etha_p_landing, CL_landing_max, T_static, theta_app, CL_landing_td) 
    return(CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si, t_cruiseh, ROC_SL_si, h_absolute_si, R_glide, S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO)
"""

"""
MTOW,OEW,Mfuelmax,PL,S,V_acc

R_eff, R_spec, R_Vmax, V_Cruise_eff, V_Cruise_spec, V_max, PowerReq, PowerAva, PowerCruise_spec, V_min, V_Loiter, V_array

"""