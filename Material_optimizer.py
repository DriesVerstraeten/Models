# -*- coding: utf-8 -*-
"""
Created on Thu Jun 01 09:32:13 2017

@author: 123
"""
#Here we put material properties, and by importing results from other files optimize according to weight-cost model. 
#Need to find this model though. 
#Inputs for each material: stiffness, cost per kg, density, other stuff if needed.

import numpy as np

rho = np.zeros(16)
E = np.zeros(16)
G = np.zeros(16)
mu = np.zeros(16)
Fty = np.zeros(16)
Ftu = np.zeros(16)
Fsu = np.zeros(16)
Fbru = np.zeros(16)

#2024-T3
rho[0] = 0.1 #density lb/in^3
E[0] = 10.5*1000 #tensile modulus ksi
G[0] = 4.0 * 1000 #shear modulus ksi
mu[0] = 0.33 #poisson
Fty[0] = 47. #yield tensile ksi
Ftu[0] = 64. #ultimate tensile ksi
Fsu[0] = 39. #ultimate shear ksi
Fbru[0] = 104. # ultimate bearing e/D=1.5 ksi
#add cost from somewhere

#2024-T4
rho[1] = 0.1 #density lb/in^3
E[1] = 10.5*1000 #tensile modulus ksi
G[1] = 4.0 * 1000 #shear modulus ksi
mu[1] = 0.33 #poisson
Fty[1] = 40. #yield tensile ksi
Ftu[1] = 62. #ultimate tensile ksi
Fsu[1] = 37. #ultimate shear ksi
Fbru[1] = 93. # ultimate bearing e/D=1.5 ksi

#6061-T6
rho[2] = 0.098 #density lb/in^3
E[2] = 9.9*1000 #tensile modulus ksi
G[2] = 3.8 * 1000 #shear modulus ksi
mu[2] = 0.33 #poisson
Fty[2] = 36. #yield tensile ksi
Ftu[2] = 42. #ultimate tensile ksi
Fsu[2] = 36. #ultimate shear ksi
Fbru[2] = 67. # ultimate bearing e/D=1.5 ksi

#7075-T6
rho[3] = 0.101 #density lb/in^3
E[3] = 10.3*1000 #tensile modulus ksi
G[3] = 3.9 * 1000 #shear modulus ksi
mu[3] = 0.33 #poisson
Fty[3] = 70. #yield tensile ksi
Ftu[3] = 78. #ultimate tensile ksi
Fsu[3] = 47. #ultimate shear ksi
Fbru[3] = 121. # ultimate bearing e/D=1.5 ksi

#AISI 1025
rho[4] = 0.284 #density lb/in^3
E[4] = 29.*1000 #tensile modulus ksi
G[4] = 11. * 1000 #shear modulus ksi
mu[4] = 0.32 #poisson
Fty[4] = 36. #yield tensile ksi
Ftu[4] = 55. #ultimate tensile ksi
Fsu[4] = 35. #ultimate shear ksi
Fbru[4] = 90. # ultimate bearing e/D=2 ksi

#AISI 4130 (t<0.118inch)
rho[5] = 0.283 #density lb/in^3
E[5] = 29.*1000 #tensile modulus ksi
G[5] = 11. * 1000 #shear modulus ksi
mu[5] = 0.32 #poisson
Fty[5] = 75. #yield tensile ksi
Ftu[5] = 95. #ultimate tensile ksi
Fsu[5] = 57. #ultimate shear ksi
Fbru[5] = 200. # ultimate bearing e/D=2 ksi

#AISI 4130 (t>0.118inch)
rho[6] = 0.283 #density lb/in^3
E[6] = 29.*1000 #tensile modulus ksi
G[6] = 11. * 1000 #shear modulus ksi
mu[6] = 0.32 #poisson
Fty[6] = 70. #yield tensile ksi
Ftu[6] = 90. #ultimate tensile ksi
Fsu[6] = 54. #ultimate shear ksi
Fbru[6] = 190. # ultimate bearing e/D=2 ksi

#AISI 4130 (t<0.118inch tempered and quenched)
rho[7] = 0.283 #density lb/in^3
E[7] = 29.*1000 #tensile modulus ksi
G[7] = 11. * 1000 #shear modulus ksi
mu[7] = 0.32 #poisson
Fty[7] = 100. #yield tensile ksi
Ftu[7] = 125. #ultimate tensile ksi
Fsu[7] = 75. #ultimate shear ksi
Fbru[7] = 146. # ultimate bearing e/D=1.5 ksi

#AISI 4340
rho[8] = 0.283 #density lb/in^3
E[8] = 29.*1000 #tensile modulus ksi
G[8] = 11. * 1000 #shear modulus ksi
mu[8] = 0.32 #poisson
Fty[8] = 217. #yield tensile ksi
Ftu[8] = 260. #ultimate tensile ksi
Fsu[8] = 156. #ultimate shear ksi
Fbru[8] = 347. # ultimate bearing e/D=1.5 ksi

#300M
rho[9] = 0.283 #density lb/in^3
E[9] = 29.*1000 #tensile modulus ksi
G[9] = 11. * 1000 #shear modulus ksi
mu[9] = 0.32 #poisson
Fty[9] = 220. #yield tensile ksi
Ftu[9] = 270. #ultimate tensile ksi
Fsu[9] = 162. #ultimate shear ksi
Fbru[9] = 414. # ultimate bearing e/D=1.5 ksi

#E-glass
rho[10] = 0.094 #density lb/in^3
E[10] = 10.*1000 #tensile modulus ksi
G[10] = 4. * 1000 #shear modulus ksi
mu[10] = 0.2 #poisson
Fty[10] = 27. #yield tensile ksi
Ftu[10] = 27. #ultimate tensile ksi
Fsu[10] = 0. #ultimate shear ksi
Fbru[10] = Ftu[10] # ultimate bearing e/D=1.5 ksi

#S-glass
rho[11] = 0.092 #density lb/in^3
E[11] = 7. * 1000 #tensile modulus ksi
G[11] = 0.6 * 1000 #shear modulus ksi
mu[11] = 0.26 #poisson
Fty[11] = 90. #yield tensile ksi
Ftu[11] = 90. #ultimate tensile ksi
Fsu[11] = 0. #ultimate shear ksi
Fbru[11] = Ftu[11] # ultimate bearing e/D=1.5 ksi

#High-modulus carbon
rho[12] = 0.072 #density lb/in^3
E[12] = 53.*1000 #tensile modulus ksi
G[12] = 2.7 * 1000 #shear modulus ksi
mu[12] = 0.2 #poisson
Fty[12] = 190. #yield tensile ksi
Ftu[12] = 190. #ultimate tensile ksi
Fsu[12] = 0. #ultimate shear ksi
Fbru[12] = Ftu[12] # ultimate bearing e/D=1.5 ksi

#High-modulus carbon
rho[13] = 0.072 #density lb/in^3
E[13] = 53.*1000 #tensile modulus ksi
G[13] = 2.7 * 1000 #shear modulus ksi
mu[13] = 0.2 #poisson
Fty[13] = 190. #yield tensile ksi
Ftu[13] = 190. #ultimate tensile ksi
Fsu[13] = 0. #ultimate shear ksi
Fbru[13] = Ftu[13] # ultimate bearing e/D=1.5 ksi

#High-strength carbon
rho[14] = 0.065 #density lb/in^3
E[14] = 35.*1000 #tensile modulus ksi
G[14] = 3.6 * 1000 #shear modulus ksi
mu[14] = 0.3 #poisson
Fty[14] = 320. #yield tensile ksi
Ftu[14] = 320. #ultimate tensile ksi
Fsu[14] = 0. #ultimate shear ksi
Fbru[14] = Ftu[14] # ultimate bearing e/D=1.5 ksi

#Aramid
rho[15] = 0.052 #density lb/in^3
E[15] = 18.*1000 #tensile modulus ksi
G[15] = 4. * 1000 #shear modulus ksi
mu[15] = 0.36 #poisson
Fty[15] = 40. #yield tensile ksi
Ftu[15] = 40. #ultimate tensile ksi
Fsu[15] = 0. #ultimate shear ksi
Fbru[15] = Ftu[15] # ultimate bearing e/D=1.5 ksi

#Convertion to normal units, lb/in^3 to kg/m^3 and ksi to Pa
rho = rho * 27.7334934 * 1000
E = E * 145.03773773 * 1000000
G = G * 145.03773773 * 1000000
Fty = Fty * 145.03773773 * 1000000
Ftu = Ftu * 145.03773773 * 1000000    
Fsu = Fsu * 145.03773773 * 1000000    
Fbru = Fbru * 145.03773773 * 1000000    
    
    