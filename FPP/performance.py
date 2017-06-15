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

