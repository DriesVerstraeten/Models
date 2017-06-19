
++++++++++++++++++++++++++++++++++++++++++++++ INPUT +++++++++++++++++++++++++++++++++++++++++++++++

This program returns the loading diagram, the state space responses and sizes the rudder and ailerons. 
Note that the functions are tested by the input parameters of the Cessna Citation. 
To run this program for the UPA, simply remove the comments (#) and use the corresponding value from the parameters.py file (p. ..)

----------------------------------------------General----------------------------------------------

Wing: 			wfn_surface_area, MAC, wing_span, wing_taper, Wing_sweep, chord_root, 
Tail: 			hor_surface_aera, ver_surface_aera, _htail_span, _vtail_span, _root_vtail, tail_arm, 
inertia: 		Ixx, Iyy, Izz, Ixz 
Weights: 		OEW_weight, OEW_loc, weight_fuel, LE_MAC, arm_fuel_tank (w.r.t nose),
Stability deriv:	f_der_long, f_der_lat, f_der_norm, m_der_pitch, m_der_yaw, m_der_roll, payl_loc, req_rolrate 

----------------------------------------------PAYLOAD----------------------------------------------

Load_c = Payload(("pilot_1", 86, 4),  --->  "payload type cruise", weight [kg], index     
                ("pilot_2", 86, 3),         
                ("pass_1",  86, 2),          
                ("pass_2",  86, 1),         
                ("bagage",  88, 0))

Load_a = = Payload(("pilot_1", 86, 4),  --->  "payload type aerobatic", weight [kg], index     
                ("pilot_2", 86, 3),         
                ("pass_1",  0, 2),          
                ("pass_2",  0, 1),         
                ("bagage",  88, 0))
                              
payl_loc_p = {1: inch(131.), 2: inch(131.),  ---> location (w.r.t nose) of the payload corresponding with the index
              3: inch(214.), 4: inch(214.),
              0: inch(338.)}

OEW_mass = OEW_cg(("wing",  692.86,   7.42), 	---> Mass distribution to determine OEW & Location
                ("engine",  692.86,   7.42), 
                ("fuselage", 692.86,  7.42), 
                ("LG_front", 692.86,   7.42), 
                ("lG_rear",  692.86,   7.42),
                ("Avionics", 692.86,   7.42))

++++++++++++++++++++++++++++++++++++++++++++++ FUNCTION ++++++++++++++++++++++++++++++++++++++++++++++

----------------------------------------------ARGUMENTS----------------------------------------------

InitCondition( _aircraft, _alti, _airspeed, _fuel_used, _angl_roll = 0, _angl_attack = 0, _angl_flpath = 0, _angl_pitch = 0, _angl_yaw = 0, _rate_roll = 0, _rate_pitch = 0, _rate_yaw = 0)

loadingdia(self, fuel_used, aero=0)			

	-fuel_used = typacally zero in this case
	aero=0, if tune to one: the loading diagram of the aerodynamic payload configuration get loaded
	+ returns: diagram and most forward and aft location

aileronSiz(self, cca_ratio, init, max_def, cd0, cla, b_outer_loc_rat)
	
	-cca_ratio = ratio aileron/wing chord
	-init =  initial conditions
	-max_def = maximum deflection: 20 degrees or 0.3490 rad
	-cd0 = increase in drag (airfoil/wing) caused by one degree of aileron deflection (should be determined by aerodynamics)
	-cla = increase in lift of the airfoil/wing caused by one degree of aileron deflection
	-b_outer_loc_rat = factor of outer aileron position of half wing span (hence 1 corresponds with the tip location)
	+ returns: aileron surface aera, outer and inner location w.r.t the center wing 

rudderSiz(self, cca_ratio, init, max_def, wake_fac, tail_eff,CL_av)

	-wake_def = factor of how much the vertical tail surface aera is in the wake during spin 
	-tail_eff = efficiency (typically 0.96) 
	-CL_av = increase in lift (along the y axis) of the vertical tail caused by one degree (slideslip) angle 
        + returns: ratio of vertical rudder span and vertical tail span (b_r/b_v), check if < 1 

responsePlot(self, init, motion, t_end, indic, rud = 0)

	- Motion: ("phugoid", "shortper", "dutch", "aperiodic"), 
	- t_end = Length time-array
	- indic: for rud = 0: ("x_vel", "aoa", "pitch", "prate")
		 for rud = 1: ("sslipa", "rolla", "rrate", "yrate")
	+ returns: reponse plots and eigenmotions/values and eigenvalues of the motions 

def Longstabfixed(self, init, cc_ratio, max_def)

----------------------------------------------RUNNING----------------------------------------------

run = Flight(_aircraft, True) 
print run_c.function(arguments)



