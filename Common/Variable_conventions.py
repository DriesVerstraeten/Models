MTOW [kg] = Max take-off weight
OEW [kg] = Operative empty weight
Mfuelmax [kg] = Maximum fuel tank capacity
PL [kg] = Payload
S [m2] = Wing surface area
A [-] = Aspect ratio
e [-] = Oswald efficiency factor
CD_0 [-] = Zero lift drag coefficient
c [m] = Mean aerodynamic chord
etha_p [-] = Propellor efficiency
etha_p_landing [-] = Propellor efficiency in landing config
x_ac [-] = Location aerodynamic center with respect to MAC
downwash [-] = depsilon/dalfa
Cmac [-] = Constant moment around aerodynamic center
CL_Ah_land [-] = Lift coeficcient aircraft-less-tail in landing config
CL_alpha_Ah_land [-] = dCL/dalpha aircraft-less-tail in landing config
CL_h_land [-] = Lift coeficcient tail in landing config
CL_alpha_h_land [-] = dCL/dalpha tail in landing config
VhV2 [-] = (Vh/V)^2
P_TO [W] = Power at  take-off



"""
def main_updown(W_TO, S, A, e, CD_0, etha_p, etha_p_landing, P_TO, P, T_static, LD_max, hcruise, CL_0, CL_alpha, alpha_TO, CL_max_TO, CL_landing_max, CL_landing_td, D_prop, D_spinner, V_c, V_H):
    CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si = Takeoff(W_TO, S, A, e, CD_0, etha_p, P_TO, CL_0, CL_alpha, alpha_TO, CL_max_TO, D_prop, D_spinner, V_c, V_H)
    t_cruiseh, ROC_SL_si, h_absolute_si = Climb(W_TO, S, A, e, CD_0, etha_p, P, LD_max, hcruise)
    R_glide = Descent(W_TO, S, A, e, CD_0, LD_max, hcruise)
    S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO = Landing(W_TO, S, A, e, CD_0, etha_p_landing, CL_landing_max, T_static, theta_app, CL_landing_td) 
    return(CL_TO, CD_TO, V_S1_si, V_LOF_si, S_TO_si, t_cruiseh, ROC_SL_si, h_absolute_si, R_glide, S_landing_si, S_ground_roll_si, S_landing_reverse_si, S_ground_roll_reverse_si, V_REF, V_SO)
"""