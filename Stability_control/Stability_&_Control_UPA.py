from math import *
import numpy as np
import control as control
import scipy.io as scio
import matplotlib.pyplot as plt 
#import Init_Parameters as p
import os
import time as tm

grav = 9.81
gas_const = 287.0
heat_cap_ratio = 1.4

def convert(val, fact, inv):
    return (fact*val) if (not inv) else (val/fact)
def inch(val, inv = False):
    fact = 0.0254
    return convert(val, fact, inv)
def pound(val, inv = False):
    fact = 0.45359237
    return convert(val, fact, inv)
def foot(val, inv = False):
    fact = 0.3048
    return convert(val, fact, inv)
def kts(val, inv = False):
    fact = 0.514444444
    return convert(val, fact, inv)
def kelvin(val, inv = False):
    return (val + 273.15) if (not inv) else (val - 273.15)

############################################################################### DEFINING FUNCTIONS ###############################################################################

class Atmosphere:
    def __init__(self): # method: class constructor
        self.isa_boundary = self.isa_boundary_fill()
    def isa_boundary_fill(self):
        boundary_alti = [0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0]
        boundary_grad = [-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002]
        boundary_temp = [288.15]
        boundary_pres = [101325.0]
        boundary_dens = [1.225]
        for num_layer in np.arange(len(boundary_alti) - 1):
            floor = {"alti": boundary_alti[num_layer], "grad": boundary_grad[num_layer], "temp": boundary_temp[num_layer], "pres": boundary_pres[num_layer], "dens": boundary_dens[num_layer]}
            floor_val = self.isa_floor(boundary_alti[num_layer + 1], floor)
            boundary_temp.append(floor_val["temp"])
            boundary_pres.append(floor_val["pres"])
            boundary_dens.append(floor_val["dens"])
        return {"alti": np.asarray(boundary_alti), "grad": np.asarray(boundary_grad), "temp": np.asarray(boundary_temp), "pres": np.asarray(boundary_pres), "dens": np.asarray(boundary_dens)}
    def isa_floor(self, alti, floor):
        alti_floor, grad, temp_floor, pres_floor, dens_floor, R = floor["alti"], floor["grad"], floor["temp"], floor["pres"], floor["dens"], gas_const
        if (grad != 0):
            temp = temp_floor + grad*(alti - alti_floor)
            pres = pres_floor*(temp/temp_floor)**(-(grav/(grad*R)))
            dens = dens_floor*(temp/temp_floor)**(-(grav/(grad*R)) - 1)
            return {"temp": temp, "pres": pres, "dens": dens}
        else:
            temp = temp_floor
            fac_exp = e**(-(grav*(alti - alti_floor)/(R*temp)))
            pres = pres_floor*fac_exp
            dens = dens_floor*fac_exp
            return {"temp": temp, "pres": pres, "dens": dens}
    def isa(self, alti):
        num_layer = np.where(self.isa_boundary["alti"] <= alti)[0][-1]
        floor = {"alti": self.isa_boundary["alti"][num_layer], "grad": self.isa_boundary["grad"][num_layer], "temp": self.isa_boundary["temp"][num_layer], "pres": self.isa_boundary["pres"][num_layer], "dens": self.isa_boundary["dens"][num_layer]}
        val = self.isa_floor(alti, floor)
        return {"temp": val["temp"], "pres": val["pres"], "dens": val["dens"]}
    def isa_temp(self, alti):
        return self.isa(alti)["temp"]
    def isa_pres(self, alti):
        return self.isa(alti)["pres"]
    def isa_dens(self, alti):
        return self.isa(alti)["dens"]

class Payload:
    def __init__(self, *args):
        self.weight, self.loc = dict(), dict()
        for payload in args:
            self.weight[payload[0]] = payload[1]
            self.loc[payload[0]] = payload[2]

class OEW_cg: 
    def __init__(self, *args): 
        self.mass, self.loc = dict(), dict()
        for OEW_cg in args: 
            self.mass[OEW_cg[0]] = OEW_cg[1]
            self.loc[OEW_cg[0]] = OEW_cg[2]

#----------------------------------------------------------------------------- Aircraft Initialisation -----------------------------------------------------------------------------#

class Aircraft:
    def __init__(self,                      # class initial state; automatically invoked during class instantiation; arguments are taken from __init__() method and stored in similarly named class attributes
                 _wfn_surface_area, _wing_chord, _wing_span, _wing_taper, _wing_sweep, _chord_root,  
                 _hor_surface_aera, _ver_surface_aera, _htail_span, _vtail_span, _root_vtail, _tail_arm, 
                 _inertia, _v_stall, _v_cruise, _OEW_mass,
                 _payload_c, _payload_a ,_weight_bempty, _weight_bempty_loc, _weight_fuel, _LE_MAC, _arm_fuel_tank,
                 _f_der_long, _f_der_lat, _f_der_norm, _m_der_pitch, _m_der_yaw, _m_der_roll, _payl_loc, _req_rolrate,  
                 _L_alpha_w, _L_alpha_ht, _m_pitch_alpha_ac, _downwash, _tailvel_eff, _ih, _stab_marg):
        self.wing = {'S': _wfn_surface_area, 'MAC': _wing_chord, 'b': _wing_span, 'taper': _wing_taper, 'sweep': _wing_sweep, 'C_r': _chord_root, 'AC_per': 0.25}     # self.wing = {'S': <wing-fuselage-nacelle combi. surface area> [m2], 'MAC': <m.a.c. chord> [m], 'b': <wingspan> [m]}
        self.tail = {'S_h': _hor_surface_aera, 'S_v': _ver_surface_aera,'b_h': _htail_span, 'b_v':_vtail_span, 'Ç_r_v': _root_vtail ,'arm': _tail_arm }
        self.coeff = {'L_alpha_w': _L_alpha_w, 'L_alpha_ht': _L_alpha_ht, 'pitch_alpha_ac': _m_pitch_alpha_ac, 'downwash': _downwash, 'tailvel_eff': _tailvel_eff, 'ih': _ih, 'SM': _stab_marg}
        self.inertia, self.v_stall, self.v_cruise = _inertia, _v_stall, _v_cruise                                                                                   # self.inertia = {"xx": <I_{xx}>, "yy": <I_{yy}>, "zz": <I_{zz}>, "xz": <I_{xz}>}
        self.payload_c, self.payload_a = _payload_c, _payload_a                                                     # self.payload = class: .weight; .loc
        self.weight_bempty, self.LE_MAC, self.weight_bempty_loc, self.weight_fuel, self.arm_fuel_tank = _weight_bempty*grav, _LE_MAC, _weight_bempty_loc, _weight_fuel*grav, _arm_fuel_tank   # self.weight_ = <weight> [N]
        self.payl_loc = _payl_loc
        self.OEW_mass_dis = _OEW_mass
        self.OEW = self.OEW_calc()
        self.f_der_long, self.f_der_lat, self.f_der_norm, self.m_der_pitch, self.m_der_yaw, self.m_der_roll = _f_der_long, _f_der_lat, _f_der_norm, _m_der_pitch, _m_der_yaw, _m_der_roll,
        self.AR = (self.wing["b"]**2)/self.wing["S"]
    def OEW_calc(self):
        segm_loc, segm_weight, OEW_moment = [], [], 0
        for segm in self.OEW_mass_dis.loc.keys(): 
            segm_loc.append(self.OEW_mass_dis.loc[segm])
        for segm in self.OEW_mass_dis.mass.keys(): 
            segm_weight.append(self.OEW_mass_dis.mass[segm]*grav)
        for i in range(len(segm_loc)): 
            OEW_moment += segm_loc[i]*segm_weight[i]
        OEW_weight = sum(segm_weight)
        OEW_cg_loc = OEW_moment/OEW_weight
        return {'OEW_weight': OEW_weight, 'OEW_cg_loc': OEW_cg_loc}
    def getWeightDry(self, aero=0): 
        payl_weight, payl_moment = 0, 0
        payl_weight_list, payl_arm_list = [],[]
        if aero == True: 
            for payl in self.payload_a.weight.keys():
                payl_weight += self.payload_a.weight[payl]*grav
                payl_moment += (self.payload_a.weight[payl]*grav)*(self.payl_loc[self.payload_a.loc[payl]])
                payl_weight_list.append(payl_weight)
                payl_arm_list.append(self.payl_loc[self.payload_a.loc[payl]])
        else: 
            for payl in self.payload_c.weight.keys():
                payl_weight += self.payload_c.weight[payl]*grav
                payl_moment += (self.payload_c.weight[payl]*grav)*(self.payl_loc[self.payload_c.loc[payl]])
                payl_weight_list.append(payl_weight)
                payl_arm_list.append(self.payl_loc[self.payload_c.loc[payl]])
        payl_arm = payl_moment/payl_weight
        weight_dry = payl_weight + self.weight_bempty
        weight_dry_moment = payl_moment + self.weight_bempty*self.weight_bempty_loc
        weight_dry_loc = weight_dry_moment/weight_dry
        return {"weight": weight_dry, "moment": weight_dry_moment, "arm": weight_dry_loc,"arm_lst": payl_arm_list, "weight_lst": payl_weight_list}

#----------------------------------------------------------------------------- S&C Functions -----------------------------------------------------------------------------#

class Flight:
    def __init__(self, _aircraft, _printActivity = False):
        self.atmos = Atmosphere()   # class instantiation: creates instance <self.atmos> of class Atmosphere
        self.aircraft = _aircraft
    def loadingdia(self, fuel_used, aero=0):   #fuel > 0: bagage > 1: pas_2 > 2: pas_1 > 3: pilot_2 > 4: pilot_1
        moment_emp, weight_emp = self.aircraft.weight_bempty*self.aircraft.weight_bempty_loc, self.aircraft.weight_bempty
        moment_fl, weight_fl = self.aircraft.arm_fuel_tank*((self.aircraft.weight_fuel/grav) - fuel_used)*grav, ((self.aircraft.weight_fuel/grav) - fuel_used)*grav
        cg_LEMAC_ld = np.array([self.aircraft.weight_bempty_loc, (moment_emp+moment_fl)/(weight_fl+weight_emp)])
        weight_LEMAC_ld = np.array([self.aircraft.weight_bempty,(weight_fl+weight_emp)])
        for i in list(self.aircraft.payl_loc):
            weight_LEMAC_ld = np.append(weight_LEMAC_ld, (self.aircraft.getWeightDry(aero)['weight_lst'][i] + weight_LEMAC_ld[i+1]))
            cg_LEMAC_ld = np.append(cg_LEMAC_ld, (self.aircraft.getWeightDry(aero)['arm_lst'][i]*self.aircraft.getWeightDry(aero)['weight_lst'][i] + cg_LEMAC_ld[i+1]*weight_LEMAC_ld[i+1])/(self.aircraft.getWeightDry(aero)['weight_lst'][i]+weight_LEMAC_ld[i+1]))
        plt.plot(cg_LEMAC_ld, weight_LEMAC_ld)
        plt.show()  
        most_forw, most_aft = min(cg_LEMAC_ld), max(cg_LEMAC_ld)
        X_cg_Fl = cg_LEMAC_ld[-1]
        return {'most_forw_ld': most_forw, 'most_aft_ld': most_aft, 'X_cg_Fl': X_cg_Fl}
    def tau(self, cca_ratio): 
        tau = -4.66*(cca_ratio**4) + 8.79*(cca_ratio**3) - 6.44*(cca_ratio**2) + 2.85*cca_ratio + 0.0316
        return tau
    def aileronSiz(self, cca_ratio, init, max_def, cd0, cla, b_outer_loc_rat):
        ail_eff = self.tau(cc_ratio) 
        M = init.airspeed/(sqrt(heat_cap_ratio*gas_const*init.isa["temp"]))
        PG, diff = (1-(M**2))**0.5, 0.1
        b_inner_loc, b_outer_loc = (b_outer_loc_rat - 0.4)*(self.aircraft.wing['b']/2), b_outer_loc_rat*(self.aircraft.wing['b']/2)              
        while True:           
            Cl_da = (cla*ail_eff*self.aircraft.wing['C_r']/(self.aircraft.wing['S']*self.aircraft.wing['b']))*((b_outer_loc**2)-(b_inner_loc**2) + (4*(self.aircraft.wing['taper']-1)/(3*self.aircraft.wing['b']))*((b_outer_loc**3)-(b_inner_loc**3)))
            cla_datcom = (2*pi*self.aircraft.AR)/(2 + sqrt(((self.aircraft.AR*PG)**2)*(1 + ((tan(self.aircraft.wing['sweep'])/PG)**2)) + 4))
            Clp = -((cla_datcom+cd0)*self.aircraft.wing['C_r']*self.aircraft.wing['b']*(1+3*self.aircraft.wing['taper']))/(24*self.aircraft.wing['S'])     
            roll_rate = -(Cl_da/Clp)*0.75*max_def*((2*init.airspeed)/self.aircraft.wing['b'])
            if (roll_rate - self.aircraft.req_rolrate) > 0: 
                b_inner_loc = b_inner_loc + 0.0001
            else: 
                b_inner_loc = b_inner_loc - 0.0001
            diff = self.aircraft.req_rolrate -  roll_rate            
            if diff < 0:
                break
        ail_surface =  (b_outer_loc - b_inner_loc)*(cca_ratio*(self.aircraft.wing['C_r']*self.aircraft.wing['taper']))
        return {'ail_aera': ail_surface, 'outer_loc': b_outer_loc, 'inner_loc': b_inner_loc} 
    def rudderSiz(self, ccr_ratio, init, max_def, wake_fac, eff,CL_av, aero=0):
        aoa_spin = init.angl_attack
        spin_rate = 0.4               # CS23 spin recovery requirement 
        inertia_trans_w = np.matrix([[(cos(aoa_spin))**2, (sin(aoa_spin))**2, -sin(2*aoa_spin)], 
                                   [(sin(aoa_spin))**2, (cos(aoa_spin))**2, -sin(2*aoa_spin)], 
                                   [0.5*sin(2*aoa_spin), -0.5*sin(2*aoa_spin), cos(2*aoa_spin)]])
        inertia_b = np.matrix([[(self.aircraft.inertia['xx']**2)*self.aircraft.getWeightDry(aero)['weight']*self.aircraft.wing['b']], 
                               [(self.aircraft.inertia['zz']**2)*self.aircraft.getWeightDry(aero)['weight']*self.aircraft.wing['b']], 
                               [(self.aircraft.inertia['xz']**2)*self.aircraft.getWeightDry(aero)['weight']*self.aircraft.wing['b']]])
        inertia_w = inertia_trans_w*inertia_b
        N_SR = ((inertia_w[0,0]*inertia_w[1,0] - (inertia_w[2,0]**2))/inertia_w[0,0])*spin_rate
        b_eR_v = (2*N_SR)/(max_def*init.isa["dens"]*(init.airspeed**2)*CL_av*self.aircraft.tail["arm"]*(1-wake_fac)*self.aircraft.tail['S_v']*eff*rud_eff)
        return {'b_ev/b_v, should be < 1': b_eR_v}   #must be < 1
    def Longstabfixed(self, init, cc_ratio, max_def): 
        X_cg_flight = self.loadingdia(init.fuel_used)["X_cg_Fl"]
        elv_rat = self.tau(cc_ratio)
        X_AC_w = self.aircraft.LE_MAC + (self.aircraft.wing['MAC']*self.aircraft.wing['AC_per'])
        V_HT = (self.aircraft.tail['S_h']*self.aircraft.tail['arm'])/(self.aircraft.wing['S']*self.aircraft.wing['MAC'])
        X_nfix = (self.aircraft.coeff['L_alpha_ht']/self.aircraft.coeff['L_alpha_w'])*(1-self.aircraft.coeff['downwash'])*self.aircraft.coeff['tailvel_eff']*((self.aircraft.tail['S_h']*self.aircraft.tail['arm'])/(self.aircraft.wing['S'])) + X_AC_w
        CL_delta_elv = (self.aircraft.tail['S_h']/self.aircraft.wing['S'])*self.aircraft.coeff['tailvel_eff']*self.aircraft.coeff['L_alpha_ht']*elv_rat
        CM_delta_elv = (self.aircraft.coeff['tailvel_eff']*(self.aircraft.tail['S_h']/self.aircraft.wing['S'])*(X_cg_flight - X_AC_w) - V_HT)*self.aircraft.coeff['L_alpha_ht']*elv_rat
        Cz_q = -2.*self.aircraft.coeff['L_alpha_ht']*self.aircraft.coeff['tailvel_eff']*V_HT
        Cm_q = -1.1*self.aircraft.coeff['L_alpha_ht']*self.aircraft.coeff['tailvel_eff']*V_HT*(self.aircraft.tail['arm']/self.aircraft.wing['MAC'])

        a = np.arange(-np.pi/18., np.pi/9., 0.001)
        ah = a*(1-self.aircraft.coeff['downwash']) + self.aircraft.coeff['ih']        
        Vair = np.arange(60.,150., 0.01)
        x =  np.arange(6.5, 8., 0.01)

        X_loc = np.array([self.loadingdia(self.aircraft.weight_fuel/grav)["most_forw_ld"], self.loadingdia(0)['most_aft_ld'], self.loadingdia((self.aircraft.weight_fuel/grav), 1)["most_forw_ld"], self.loadingdia(0,1)['most_aft_ld']])
        deff = np.array([max_def,0., -max_def])

        Cm_a = []
        delta_ev_v = []
        
        for i in X_loc: 
            for j in deff: 
                Cm_delta_i = self.aircraft.coeff['pitch_alpha_ac'] +  self.aircraft.coeff['L_alpha_w']*a*((i-X_AC_w)/self.aircraft.wing['MAC']) - (self.aircraft.coeff['L_alpha_ht']*ah - CL_delta_elv*j)*(self.aircraft.coeff['tailvel_eff']*V_HT)
                Cm_a.append(Cm_delta_i)
            delta_ev = (-1/CM_delta_elv)*(self.aircraft.coeff['pitch_alpha_ac'] + (init.weight/(0.5*init.isa["dens"]*(Vair**2)*self.aircraft.wing['S']))*((i - X_nfix)/self.aircraft.wing["MAC"]))
            delta_ev_v.append(delta_ev)

        delta_elv_dn_pull_0 = (-1/CM_delta_elv)*(init.weight/(0.5*self.atmos.isa_dens(0.)*(init.airspeed**2)*self.aircraft.wing['S']))*((Cm_q/((2*init.weight)/(grav*self.atmos.isa_dens(0.)*self.aircraft.wing['MAC']*self.aircraft.wing['S']))) + ((x - X_nfix)/self.aircraft.wing['MAC']))
        delta_elv_dn_pull_c = (-1/CM_delta_elv)*(init.weight/(0.5*init.isa['dens']*(init.airspeed**2)*self.aircraft.wing['S']))*((Cm_q/(2*init.mass_coef['MAC'])) + ((x - X_nfix)/self.aircraft.wing['MAC']))

        X_pos = np.append(X_loc, X_nfix)
        labels1 = ['CG range (cruise)', 'CG range (aero)', 'Neutral point']
        points = ['d','s','o']
        plt.subplot(221)
        plt.plot((x - self.aircraft.LE_MAC)/self.aircraft.wing['MAC'], (delta_elv_dn_pull_0/np.pi)*180, label='H= 0 ft')
        plt.plot((x - self.aircraft.LE_MAC)/self.aircraft.wing['MAC'], (delta_elv_dn_pull_c/np.pi)*180, label='H= 18000 ft')        
        plt.plot([(X_pos[0] - self.aircraft.LE_MAC)/self.aircraft.wing['MAC'], (X_pos[1] - self.aircraft.LE_MAC)/self.aircraft.wing['MAC']], [0.05,0.05], marker=points[0], label=labels1[0])
        plt.plot([(X_pos[2] - self.aircraft.LE_MAC)/self.aircraft.wing['MAC'], (X_pos[3] - self.aircraft.LE_MAC)/self.aircraft.wing['MAC']], [-0.05,-0.05], marker=points[1], label=labels1[1])
        plt.plot((X_pos[4] - self.aircraft.LE_MAC)/self.aircraft.wing['MAC'], [0.], marker=points[2], label=labels1[2])
        plt.gca().invert_yaxis()
        plt.xlabel('$x_{cg}/MAC$ [-]')
        plt.ylabel('$d{\delta_e}/dn$ $[{\degree}]$')
        plt.grid(True, which='both')    
        plt.legend(loc='upper right')

        plt.subplot(222)
        labels2 = ['CG Front (cruise)', 'CG Aft (cruise & aero)', 'CG Front (aero)', 'CG Aft (aero)']
        for i in range(len(delta_ev_v)-1):
            plt.plot(Vair, delta_ev_v[i], label=labels2[i])
        plt.xlabel('V [m/s]')
        plt.ylabel('${\delta_e}$ $[{\degree}]$')
        plt.grid(True, which='both')
        plt.legend(loc='lower right') 

        plt.subplot(223)
        labels3 = ['CG Front (cruise)', 'CG Aft (cruise & aero)', 'CG Front (aero)', 'CG Aft (aero)']
        colors = ['b', 'g', 'r']
        for i in range(len(Cm_a)/4):
            plt.plot((a/np.pi)*180, Cm_a[1 + 3*i], label=labels3[i])
            plt.fill_between((a/np.pi)*180, Cm_a[1 + 3*i] + abs(Cm_a[0]-Cm_a[1]), Cm_a[1 + 3*i] - abs(Cm_a[0]-Cm_a[1]), color=colors[i], alpha=0.3)
        plt.grid(True, which='both')
        plt.xlabel('alpha [${\degree}$]')
        plt.ylabel('$C_m$ [-]')
        plt.legend(loc='lower left') 

        plt.show() 
        return{'CL_delta_elv': CL_delta_elv, "CM_delta_elv": CM_delta_elv, 'X_nfix': X_nfix }
    def stablecheck(self, fuel_used, cc_ratio, max_def):
        if (self.Longstabfixed(init, cc_ratio, max_def)['X_nfix'] > self.loadingdia(0)['most_aft_ld']) or (self.Longstabfixed(init, cc_ratio, max_def)['X_nfix'] > self.loadingdia(self.aircraft.weight_fuel)['most_aft_ld']): 
            return 'Aircraft is unstable!'
        else: 
            return 'Aircraft is stable!'
    def space_state_sys(self, init, motion, t_end, rud = 0):      #motion: "shortper", "phugoid", "aperiotic", "dutch", t_end = end time of array
        sym_dict = {"shortper": True, "phugoid": True, "dutch": False, "aperiodic": False}
        sym = sym_dict[motion]
        mass_coeff_c, mass_coeff_b, C_L, C_z_0, C_x_0 = init.mass_coef['MAC'], init.mass_coef["b"], init.zero_coef["L"], init.zero_coef["z"], init.zero_coef["x"]
        V = init.airspeed
        T = np.linspace(0, t_end , num = (t_end/0.001)+1.) 
        if sym: 
            C1 = np.matrix([[-2.*mass_coeff_c*(self.aircraft.wing['MAC']/(V**2)), 0 , 0, 0], 
                            [0, (self.aircraft.f_der_norm["alpha_dot"] - 2.*mass_coeff_c)*(self.aircraft.wing['MAC']/V), 0, 0],
                            [0, 0, (-1*self.aircraft.wing['MAC']/V), 0], 
                            [0, self.aircraft.m_der_pitch["alpha_dot"]*(self.aircraft.wing['MAC']/V), 0, -2.*mass_coeff_c*(self.aircraft.inertia["yy"]**2.)*((self.aircraft.wing['MAC']/V)**2.)]])
            C2 = np.matrix([[self.aircraft.f_der_long["u"]/V, self.aircraft.f_der_long["alpha"], C_z_0, self.aircraft.f_der_long["q"]*(self.aircraft.wing['MAC']/V)], 
                            [self.aircraft.f_der_norm["u"]/V, self.aircraft.f_der_norm["alpha"], -C_x_0, (self.aircraft.f_der_norm["q"] + 2*mass_coeff_c)*(self.aircraft.wing['MAC']/V)], 
                            [0, 0, 0, (self.aircraft.wing['MAC']/V)],
                            [self.aircraft.m_der_pitch["u"]/V, self.aircraft.m_der_pitch["alpha"], 0, self.aircraft.m_der_pitch["q"]*(self.aircraft.wing['MAC']/V)]])
            C3 = np.matrix([[self.aircraft.f_der_long["delta"]], 
                            [self.aircraft.f_der_norm["delta"]], 
                            [0], 
                            [self.aircraft.m_der_pitch["delta"]]])
            X0 = init.symm_vect
        else: 
            C1 = np.matrix([[(self.aircraft.f_der_lat["beta_dot"] - 2*mass_coeff_b)*(self.aircraft.wing["b"]/V), 0, 0, 0], 
                            [0, -1*(self.aircraft.wing["b"]/(2*V)), 0, 0,], 
                            [0, 0, -2*mass_coeff_b*(self.aircraft.inertia["xx"]**2)*((self.aircraft.wing["b"]/V)**2), 2*mass_coeff_b*(self.aircraft.inertia["xz"])*((self.aircraft.wing["b"]/V)**2)], 
                            [self.aircraft.m_der_yaw["beta_dot"]*(self.aircraft.wing["b"]/V), 0, 2*mass_coeff_b*(self.aircraft.inertia["xz"])*((self.aircraft.wing["b"]/V)**2), -2*mass_coeff_b*(self.aircraft.inertia["zz"]**2)*((self.aircraft.wing["b"]/V)**2)]])
            
            C2 = np.matrix([[self.aircraft.f_der_lat["beta"], C_L, self.aircraft.f_der_lat["p"]*(self.aircraft.wing["b"]/(2*V)), (self.aircraft.f_der_lat["r"] - 4*mass_coeff_b)*(self.aircraft.wing["b"]/(2*V))], 
                            [0, 0, (self.aircraft.wing["b"]/(2*V)), 0], 
                            [self.aircraft.m_der_roll["beta"], 0, self.aircraft.m_der_roll["p"]*(self.aircraft.wing["b"]/(2*V)), self.aircraft.m_der_roll["r"]*(self.aircraft.wing["b"]/(2*V))], 
                            [self.aircraft.m_der_yaw["beta"], 0, self.aircraft.m_der_yaw["p"]*(self.aircraft.wing["b"]/(2*V)), self.aircraft.m_der_yaw["r"]*(self.aircraft.wing["b"]/(2*V))]])
            
            C3 = np.matrix([[self.aircraft.f_der_lat["delta_a"], self.aircraft.f_der_lat["delta_r"]],
                            [0, 0], 
                            [self.aircraft.m_der_roll["delta_a"], self.aircraft.m_der_roll["delta_r"]], 
                            [self.aircraft.m_der_yaw["delta_a"], self.aircraft.m_der_yaw["delta_r"]]])
            X0 = init.asym_vect
        A = np.linalg.inv(C1)*(-1.*C2)    
        B = np.linalg.inv(C1)*(-1.*C3)
        C = np.eye(A.shape[1])
        D = np.zeros(shape=(B.shape))
        U1, U2 = np.full(((1/0.001)+1),(-pi/180.)), np.full(((t_end/0.001) - (1/0.001)),0.0)  
        U = np.append(U1,U2)
        sys = control.ss(A,B,C,D)
        eigv = np.linalg.eig(A)
        if sym ==1:             
            res = control.forced_response(sys,T,U,X0)        #Returns: T-array, yout-array[x_vec], xout-array
            res[1][0] = res[1][0] + V
        elif sym == 0 and rud == 1:  
            sys = control.ss(A,B[:,1],C, np.zeros(shape=(4,1)))
            res = control.forced_response(sys,T,U,X0)
        elif sym == 0 and rud == 0:
            res = control.initial_response(sys,T,X0,rud) 
        return  res[0], res[1][0], res[1][1], res[1][2], res[1][3], eigv   # T_arr, x_vel, aoa, pitch, prate, sslipa, rolla, rrate, yrate
    def responsePlot(self, init, motion, t_end, indic, rud = 0):
        vec_indic = {"x_vel": 1, "aoa": 2, "pitch": 3, "prate": 4,
                     "sslipa": 1, "rolla": 2, "rrate": 3, "yrate": 4}
        resp = self.space_state_sys(init, motion, t_end, rud)
        plt.plot(resp[0], resp[vec_indic[indic]])
        plt.show()
        return {"eigenvalue": resp[5][0] }

#----------------------------------------------------------------------------- Initial condition functions -----------------------------------------------------------------------------#
   
class InitCondition:
    def __init__(self, _aircraft, _alti, _airspeed, _fuel_used, aero = 0, _angl_roll = 0, _angl_attack = 0, _angl_flpath = 0, _angl_pitch = 0, _angl_yaw = 0, _rate_roll = 0, _rate_pitch = 0, _rate_yaw = 0):
        self.aircraft, self.atmos = _aircraft, Atmosphere()
        self.alti, self.airspeed, self.fuel_used, self.angl_attack, self.angl_flpath, self.angl_roll, self.angl_pitch, self.angl_yaw, self.rate_roll, self.rate_pitch, self.rate_yaw = _alti, _airspeed, _fuel_used, _angl_attack, _angl_flpath, _angl_roll, _angl_pitch, _angl_yaw, _rate_roll, _rate_pitch, _rate_yaw
        self.weight = self.aircraft.getWeightDry(aero)['weight'] + self.aircraft.weight_fuel - self.fuel_used
        self.isa = self.atmos.isa(self.alti)
        self.symm_vect, self.asym_vect = self.getSymmVect(), self.getAsymVect()
        self.mass_coef, self.zero_coef = self.getMassCoef(), self.getZeroCoef()
    def getSymmVect(self):
        return np.matrix([[0.],
                            [self.angl_attack],
                            [self.angl_pitch],
                            [self.rate_pitch]])
    def getAsymVect(self):
        return np.matrix([[self.angl_yaw],
                            [self.angl_roll],
                            [self.rate_roll],
                            [self.rate_yaw]])
    def getMassCoef(self):
        return {'MAC': (self.weight/grav)/(self.isa["dens"]*self.aircraft.wing['S']*self.aircraft.wing['MAC']), 'b': (self.weight/grav)/(self.isa["dens"]*self.aircraft.wing['S']*self.aircraft.wing['b'])}
    def getZeroCoef(self): 
        return {'L': (self.weight*cos(self.angl_flpath))/(0.5*self.isa["dens"]*(self.airspeed**2.)*self.aircraft.wing["S"]), 'z': (-1*self.weight*cos(self.angl_pitch))/(0.5*self.isa["dens"]*(self.airspeed**2.)*self.aircraft.wing["S"]), 'x': (self.weight*sin(self.angl_pitch))/(0.5*self.isa["dens"]*(self.airspeed**2)*self.aircraft.wing["S"])}   


############################################################################### DEFINING PARAMETERS ###############################################################################

#-----------------------------------------------------------------------------  CESSNA CITATION 550 -----------------------------------------------------------------------------#
                 
_f_der_long_c = {'u': -0.0279, "alpha": -0.4797, "alpha_dot": 0.0833, 'q': -0.2817, "delta": -0.0373}                                 # longitudinal force derivatives, i.e., C_{X}
_f_der_norm_c = {'u': -0.3762, "alpha": -5.7434, "alpha_dot": -0.0035, 'q': -5.6629, "delta": -0.6961}                                # normal force derivatives, i.e., C_{Z}
_f_der_lat_c = {"beta": -0.75, "beta_dot": 0.0, 'p': -0.0304, 'r': 0.8495, "delta_a": -0.04, "delta_r": 0.23}                         # lateral force derivatives, i.e., C_{Y}
_m_der_pitch_c = {'0': 0.0297, 'u': 0.0699, "alpha": -0.5626, "alpha_dot": 0.178, 'q': -8.7941, "delta": -1.1642, "T_c": -0.0064}     # pitch moment derivatives, i.e., C_{m}
_m_der_yaw_c = {"beta": 0.1348, "beta_dot": 0.0, 'p': -0.0602, 'r': -0.2061, "delta_a": -0.012, "delta_r": -0.0939}                   # yaw moment derivatives, i.e., C_{n}
_m_der_roll_c = {"beta": -0.1026, 'p': -0.7108, 'r': 0.2376, "delta_a": -0.2309, "delta_r": 0.0344}                                   # roll moment derivatives, i.e., C_{l}

load_c = Payload(("pilot_1", 90, 1),          
               ("pilot_2", 90, 2),          
               ("pass_1", 86, 9),           
               ("pass_2", 79, 4),             
               ("pass_3", 87, 7),           
               ("pass_4", 90, 8),
               ("pass_5", 79, 5),
               ("pass_6", 78, 3),
               ("pass_7", 75, 6), 
               ("bagage", 50, 0))

load_a = Payload(("pilot_1", 0, 1),          
               ("pilot_2", 0, 2),          
               ("pass_1", 86, 9),           
               ("pass_2", 0, 4),             
               ("pass_3", 90, 7),           
               ("pass_4", 90, 8),
               ("pass_5", 0, 5),
               ("pass_6", 0, 3),
               ("pass_7", 75, 6), 
               ("bagage", 50, 0))

payl_loc_c = {1: inch(131.), 2: inch(131.),
              3: inch(214.), 4: inch(214.),
              5: inch(251),  6: inch(251.),
              7: inch(288.), 8: inch(288.),
              9: inch(170.), 0: inch(325.)}

OEW_mass = OEW_cg(("wing",  100.,   3.25), 
                ("engine",  150.,   1.75), 
                ("fuselage", 450.,  4.5), 
                ("LG_front", 50.,   1.75), 
                ("lG_rear",  80.,   3.75 ),
                ("Avionics", 20.,   1.75))

c550 = Aircraft(30., 2.0569, 15.911, 0.316, 0., 2.15,                 
                6.23, 1.8 , 5.79, 1.5, 1., 5,                             
                {"xx": sqrt(0.019), "yy": sqrt(1.3925), "zz": sqrt(0.042), "xz": 0.002}, 62, 92.6, OEW_mass,  # class instantiation
                load_c, load_a, pound(9165.0), inch(292.18), 1814.39, inch(261.), inch(300), 
                _f_der_long_c, _f_der_lat_c, _f_der_norm_c, _m_der_pitch_c, _m_der_yaw_c, _m_der_roll_c, payl_loc_c, 0.0175,
                4.1, 3.91, 0.0197, 0.45, 0.9, 0., 0.05)

                #_wfn_surface_area, _wing_chord, _wing_span, _wing_taper, _wing_sweep, _chord_root,  
                #_hor_surface_aera, _ver_surface_aera, _htail_span, _vtail_span, _root_vtail, _tail_arm, 
                #_inertia, _v_stall, _v_cruise, _OEW_mass,
                #_payload_c, _payload_a ,_weight_bempty, _weight_bempty_loc, _weight_fuel, _LE_MAC, _arm_fuel_tank,
                #_f_der_long, _f_der_lat, _f_der_norm, _m_der_pitch, _m_der_yaw, _m_der_roll, _payl_loc, _req_rolrate,  
                #_L_alpha_w, _L_alpha_ht, _m_pitch_alpha_ac, _downwash, _tailvel_eff, _ih, _stab_marg
#----------------------------------------------------------------------------- UPA -----------------------------------------------------------------------------#

#_f_der_long_p = {'u': -0.0279, "alpha": -0.4797, "alpha_dot": 0.0833, 'q': -0.2817, "delta": -0.0373}                                 # longitudinal force derivatives, i.e., C_{X}
#_f_der_norm_p = {'u': -0.3762, "alpha": -5.7434, "alpha_dot": -0.0035, 'q': -5.6629, "delta": -0.6961}                                # normal force derivatives, i.e., C_{Z}
#_f_der_lat_p = {"beta": -0.75, "beta_dot": 0.0, 'p': -0.0304, 'r': 0.8495, "delta_a": -0.04, "delta_r": 0.23}                         # lateral force derivatives, i.e., C_{Y}
#_m_der_pitch_p = {'0': 0.0297, 'u': 0.0699, "alpha": -0.5626, "alpha_dot": 0.178, 'q': -8.7941, "delta": -1.1642, "T_c": -0.0064}     # pitch moment derivatives, i.e., C_{m}
#_m_der_yaw_p = {"beta": 0.1348, "beta_dot": 0.0, 'p': -0.0602, 'r': -0.2061, "delta_a": -0.012, "delta_r": -0.0939}                   # yaw moment derivatives, i.e., C_{n}
#_m_der_roll_p = {"beta": -0.1026, 'p': -0.7108, 'r': 0.2376, "delta_a": -0.2309, "delta_r": 0.0344}                                   # roll moment derivatives, i.e., C_{l}

#Load_p = Payload(("pilot_1", 86, 4),        # class instantiation: creates instance <aerotica> of class Crew.
#               ("pilot_2", 86, 3),          # arguments: tuple (<name>, <weight (kg)>, <payl_loc>) for each crew member
#               ("pass_1",  86, 2),          #   <payl_loc> = [1, 10]
#               ("pass_2",  86, 1),          # arguments: tuple (<name>, <weight (kg)>, <payl_loc>) for each piece of bagage           
#               ("bagage",  88, 0))
                              
#payl_loc_p = {1: 1.2, 2: 1.2,  # location of the pilots, passangers and bagage w.r.t the datum line or LE_MAC
#              3: 2.4, 4: 2.4,
#              0: 3}

#UPA = Aircraft(p.S, p.MAC, p.b, (p.c_t/p.c_r), 0., p.c_r,                 
#                p.S_ht, 1.8 , p.S_ht, 1.5, 1., 5.,                             
#               {"xx": sqrt(0.019), "yy": sqrt(1.3925), "zz": sqrt(0.042), "xz": 0.002},  # class instantiation
#                load_c, pound(9165.0), inch(292.18), pound(4000.), inch(261.), inch(300),
#                _f_der_long_c, _f_der_lat_c, _f_der_norm_c, _m_der_pitch_c, _m_der_yaw_c, _m_der_roll_c, payl_loc_c, 0.0175)

                 #_wfn_surface_area, _wing_chord, _wing_span, _hor_surface_aera, _wing_taper, _wing_sweep, _chord_root, 
                 #_hor_surface_aera, _ver_surface_aera, _htail_span, _vtail_span, _root_vtail, tail_arm, 
                 #_inertia,
                 #_payload, _weight_bempty, _weight_bempty_loc, _weight_fuel, _LE_MAC, _arm_fuel_tank,
                 #_f_der_long, _f_der_lat, _f_der_norm, _m_der_pitch, _m_der_yaw, _m_der_roll, _payl_loc, _req_rolrate    

#----------------------------------------------------------------------------- Payload -----------------------------------------------------------------------------#

############################################################################### RUNNING FUNCTIONS ###############################################################################

#---------------------------------------------------------------CESSNA CITATION---------------------------------------------------------------#

run_c = Flight(c550, True)                                          
init_test_c550 = InitCondition(c550, 3000., 93., 900., 0.)    #_aircraft, _alti, _airspeed, _fuel_used, _angl_roll = 0, _angl_attack = 0, _angl_flpath = 0, _angl_pitch = 0, _angl_yaw = 0, _rate_roll = 0, _rate_pitch = 0, _rate_yaw = 0          

init_rud_c = InitCondition(c550, 5000., 44., 600., 0., 1)     # V_stall
init_cruise_c = InitCondition(c550, 10000., 150., 600.)


#aileron_size =  run_c.aileronSiz(0.3, init_cruise_C, 0.349066, 0.010, 0.05524, 0.9)         #cca_ratio, init, max_def, cd0, cla, b_outer_loc_rat
#rud_size = run_c.rudderSiz(0.4, init_rud_c, 0.349066, 0, 0.96, 5.4 )                      #cca_ratio, init, max_def, wake_fac, tail_eff,CL_av
#loading_dia =  run_c.loadingdia(1814.369, 0)
#response_plot =  run_c.responsePlot (init_cruise_c, "phugoid",200, "x_vel", 1)    # init, motion, t_end, indic, rud = 0

#("phugoid", "shortper", "dutch", "aperiodic"), Length time-array, plot: ("x_vel", "aoa", "pitch", "prate"), 1 = rudder input, 0 = aileron input
#                                                                        ("sslipa", "rolla", "rrate", "yrate")

longstab = run_c.Longstabfixed(init_cruise_c, 0.4, 0.349066)
#---------------------------------------------------------------UPA---------------------------------------------------------------#

#run_c = Flight(UPA, True)                                                                       
#init_test_UPA = InitCondition(UPA, 3000., 93., 900., 0.)    #_aircraft, _alti, _airspeed, _fuel_used, _angl_roll = 0, _angl_attack = 0, _angl_flpath = 0, _angl_pitch = 0, _angl_yaw = 0, _rate_roll = 0, _rate_pitch = 0, _rate_yaw = 0

#init_rud_p = InitCondition(UPA, 5000., 34., 50. 0. 1.)                                     #random altitude, V_stall ,fuel used does't matter,
#init_cruise_p = InitCondition(UPA, p.h_cruise, p.V_cruise., 50.) 

#aileron_size =  run_p.aileronSiz(0.3, init_test_ail_p, 0.349066, 0.010, 0.05524, 0.9)           #cca_ratio, init, max_def, cd0, cla, b_outer_loc_rat
#rud_size = run_p.rudderSiz(0.4, init_test_rud_p, 0.349066, 0.4, 0.96, 4.4 )                     #cca_ratio, init, max_def, wake_fac, tail_eff,CL_av
#loading_dia =  run_p.loadingdia(0.)
#response_plot =  run_c.responsePlot (init_cruise, "phugoid",200, "x_vel", 1) 

#---------------------------------------------------------------Printing---------------------------------------------------------------#

print longstab

#----------------------------- Required inputs for UPA -----------------------------# 
# - Position and mass of the pilots, passengers & bagage 
# - Stability derivatives 
# - OEW with cg Location, Begin Fuel Weight and moment caused by fuel at
# - Wing area, Wing  chord, Wing span, drag coeff 0, lift coeff aoa slope,  oswald fact, inertia, datum (distance to LE_MAC)
