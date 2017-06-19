# -*- coding: utf-8 -*-

from math import *
import numpy as np
import control as control
import scipy.io as scio
import matplotlib.pyplot as plt 
import StabilityControlUPA as sc
import os
import time as tm

                        
                            
StabCheck = sc.Flight(sc.Aircraft(30., 2.0569, 15.911, 0.316, 0., 2.15,                 
                6.23, 1.8 , 5.79, 1.5, 1., 5, 1.2,                         
                {"xx": sqrt(0.019), "yy": sqrt(1.3925), "zz": sqrt(0.042), "xz": 0.002}, 62, 92.6, 

                sc.OEW_cg(("wing",  692.86,   7.42), 
                          ("engine",  692.86,   7.42), 
                          ("fuselage", 692.86,  7.42), 
                          ("LG_front", 692.86,   7.42), 
                          ("lG_rear",  692.86,   7.42),
                          ("Avionics", 692.86,   7.42)),

                sc.Payload(("pilot_1", 90, 1),          
                           ("pilot_2", 90, 2),          
                           ("pass_1", 86, 9),           
                           ("pass_2", 79, 4),             
                           ("pass_3", 87, 7),           
                           ("pass_4", 90, 8),
                           ("pass_5", 79, 5),
                           ("pass_6", 78, 3),
                           ("pass_7", 75, 6), 
                           ("bagage", 50, 0)), 

                sc.Payload(("pilot_1", 0, 1),          
                            ("pilot_2", 0, 2),          
                            ("pass_1", 86, 9),           
                            ("pass_2", 0, 4),             
                            ("pass_3", 90, 7),           
                            ("pass_4", 90, 8),
                            ("pass_5", 0, 5),
                            ("pass_6", 0, 3),
                            ("pass_7", 75, 6), 
                            ("bagage", 50, 0)), 

                1814.39, sc.inch(261.), sc.inch(300),

                {'u': -0.0279, "alpha": -0.4797, "alpha_dot": 0.0833, 'q': -0.2817, "delta": -0.0373}, 
                {'u': -0.3762, "alpha": -5.7434, "alpha_dot": -0.0035, 'q': -5.6629, "delta": -0.6961},
                {"beta": -0.75, "beta_dot": 0.0, 'p': -0.0304, 'r': 0.8495, "delta_a": -0.04, "delta_r": 0.23}, 
                {'0': 0.0297, 'u': 0.0699, "alpha": -0.5626, "alpha_dot": 0.178, 'q': -8.7941, "delta": -1.1642, "T_c": -0.0064}, 
                {"beta": 0.1348, "beta_dot": 0.0, 'p': -0.0602, 'r': -0.2061, "delta_a": -0.012, "delta_r": -0.0939}, 
                {"beta": -0.1026, 'p': -0.7108, 'r': 0.2376, "delta_a": -0.2309, "delta_r": 0.0344}, 

                {1: sc.inch(131.), 2: sc.inch(131.),
                 3: sc.inch(214.), 4: sc.inch(214.),
                 5: sc.inch(251),  6: sc.inch(251.),
                 7: sc.inch(288.), 8: sc.inch(288.),
                 9: sc.inch(170.), 0: sc.inch(325.)}, 

                0.0175, 4.1, 3.91, 0.0197, 0.45, 0.9, 0., 0.05)).Longstabfixed(sc.InitCondition(sc.Aircraft(30., 2.0569, 15.911, 0.316, 0., 2.15,                 
                                                                                                6.23, 1.8 , 5.79, 1.5, 1., 5, 1.2,                             
                                                                                                {"xx": sqrt(0.019), "yy": sqrt(1.3925), "zz": sqrt(0.042), "xz": 0.002}, 62, 92.6, 

                                                                                                sc.OEW_cg(("wing",  100.,   3.25), 
                                                                                                          ("engine",  150.,   1.75), 
                                                                                                          ("fuselage", 450.,  4.5), 
                                                                                                          ("LG_front", 50.,   1.75), 
                                                                                                          ("lG_rear",  80.,   3.75 ),
                                                                                                          ("Avionics", 20.,   1.75)),  

                                                                                                sc.Payload(("pilot_1", 90, 1),          
                                                                                                           ("pilot_2", 90, 2),          
                                                                                                           ("pass_1", 86, 9),           
                                                                                                           ("pass_2", 79, 4),             
                                                                                                           ("pass_3", 87, 7),           
                                                                                                           ("pass_4", 90, 8),
                                                                                                           ("pass_5", 79, 5),
                                                                                                           ("pass_6", 78, 3),
                                                                                                           ("pass_7", 75, 6), 
                                                                                                           ("bagage", 50, 0)), 

                                                                                                sc.Payload(("pilot_1", 90, 1),          
                                                                                                           ("pilot_2", 90, 2),          
                                                                                                           ("pass_1", 86, 9),           
                                                                                                           ("pass_2", 79, 4),             
                                                                                                           ("pass_3", 87, 7),           
                                                                                                           ("pass_4", 90, 8),
                                                                                                           ("pass_5", 79, 5),
                                                                                                           ("pass_6", 78, 3),
                                                                                                           ("pass_7", 75, 6), 
                                                                                                           ("bagage", 50, 0)), 

                                                                                                1814.39, sc.inch(261.), sc.inch(300),

                                                                                                {'u': -0.0279, "alpha": -0.4797, "alpha_dot": 0.0833, 'q': -0.2817, "delta": -0.0373}, 
                                                                                                {'u': -0.3762, "alpha": -5.7434, "alpha_dot": -0.0035, 'q': -5.6629, "delta": -0.6961},
                                                                                                {"beta": -0.75, "beta_dot": 0.0, 'p': -0.0304, 'r': 0.8495, "delta_a": -0.04, "delta_r": 0.23}, 
                                                                                                {'0': 0.0297, 'u': 0.0699, "alpha": -0.5626, "alpha_dot": 0.178, 'q': -8.7941, "delta": -1.1642, "T_c": -0.0064}, 
                                                                                                {"beta": 0.1348, "beta_dot": 0.0, 'p': -0.0602, 'r': -0.2061, "delta_a": -0.012, "delta_r": -0.0939}, 
                                                                                                {"beta": -0.1026, 'p': -0.7108, 'r': 0.2376, "delta_a": -0.2309, "delta_r": 0.0344}, 

                                                                                                {1: sc.inch(131.), 2: sc.inch(131.),
                                                                                                 3: sc.inch(214.), 4: sc.inch(214.),
                                                                                                 5: sc.inch(251),  6: sc.inch(251.),
                                                                                                 7: sc.inch(288.), 8: sc.inch(288.),
                                                                                                 9: sc.inch(170.), 0: sc.inch(325.)}, 

                                                                                                0.0175, 4.1, 3.91, 0.0197, 0.45, 0.9, 0., 0.05), 10000., 150., 600.), 0.5, 0.4363)
          
print StabCheck                                                                                                     