import math
import numpy as np

def class2est(inputlist,inputlist2,mtow,oew,fuelw,maxfuelw):

    Cd0 = inputlist[0]
    e = inputlist[1]
    A = inputlist[2]
    propefficiency = inputlist[3]
    u = inputlist[4]
    Vcruise = inputlist[5]
    height = inputlist[6]
    payload = inputlist[7]
    flightrange = inputlist[8]
    loitertime = inputlist[9]

    aref = inputlist[10]
    bref = inputlist[11]

    CLmaxLand = inputlist[12]
    Landingdistance = inputlist[13]

    Cp = u/5918352.

    Vland = np.sqrt(Landingdistance/0.5915)

    WSland = (0.5*CLmaxLand*Vland*Vland*1.225)

    S = mtow/WSland*9.81
    
    sweepquarter = inputlist2[0]    #F4
    taper = inputlist2[1]           #F5
    tc = inputlist2[2]              #F6
    Nult = inputlist2[3]            #F7
    rootcoord = inputlist2[4]       #F11
    Vh = inputlist2[5]              #F12
    N_pax = inputlist2[6]           #F13
    M_D = inputlist2[7]             #F14
    Sh = inputlist2[8]              #F16
    Ah = inputlist2[9]              #F17
    Ch_max = inputlist2[10]         #F18
    trh = inputlist2[11]            #F19
    bh = inputlist2[12]             #F20
    lh = inputlist2[13]             #F21
    Sv = inputlist2[14]             #F22
    Av = inputlist2[15]             #F23
    trv = inputlist2[16]            #F24
    bv = inputlist2[17]             #F25
    cv_max = inputlist2[18]         #F26
    fuselage_l = inputlist2[19]     #J3
    fuse_width = inputlist2[20]     #J4
    fuse_height = inputlist2[21]    #J5
    Vc = inputlist2[22]             #J6

    sep_tanks = inputlist2[23]      #J9
    We = inputlist2[24]             #J10
    Ne = inputlist2[25]             #J11
    
    Nult1 = inputlist2[26]          #J13
    lsm = inputlist2[27]            #J14
    W_L = inputlist2[28]            #J15
    integralfrac = inputlist2[29]

    W_wing = (96.948*((((mtow*2.204)*Nult/(10.0**5))**0.65)*((A/(math.cos(sweepquarter)))**0.57)*(((S/(0.3048*0.3048))/100)**0.61)*(((1+taper)/(2*tc))**0.36)*((1+(Vh/0.5144444/500.0))**0.5))**0.993)*0.453592
    W_h = (127*(((((mtow*2.204)*Nult/(10**5))**0.87)*(((Sh/(0.3048*0.3048))/100.0)**1.2)*0.289*(((lh/0.3048)/10)**0.483)*(((bh/0.3048)/(trh/0.3048))**0.5))**0.458))*0.45359237
    W_v =(2*(98.5*(((((mtow*2.204)*Nult/(10.0**5))**0.87)*(((Sv/(0.3048*0.3048))/100.0)**1.2)*0.289*(((bv/0.3048)/(trv/0.3048))**0.5))**0.458)))*0.45359237

    W_fuselage = (200*((mtow*2.204)*Nult/(10.0**5))**0.286*((fuselage_l/0.3048)/10)**0.857*(((fuse_width+fuse_height)/0.3048)/10.0*((Vc/0.51444444)/100.0)**0.338)**1.1)/2.204
    W_gear =(0.054*((lsm/.3048)**0.501)*(((mtow*2.204)*Nult1)**0.684))*0.45359237
    W_pwr =2.575*Ne*(We*2.204)**0.922/2.204
    W_Fuel_sys =(2.49*((maxfuelw*2.204/6.55)**0.6*(1/(1+integralfrac))**0.3*(sep_tanks)**0.2*(Ne)**0.13)**1.21)/2.204
    W_flightcontrol =1.066*((mtow*2.204)**0.626)*0.453592
    W_iae = (33.0*N_pax)*0.453592
    W_els = 426.0*((((W_Fuel_sys/0.45359237)+(W_iae/0.45359237))/1000)**0.51)*0.45359237
    W_api = 0.265*((mtow*2.204)**0.52)*(N_pax**0.68)*((W_iae/0.453592)**0.17)*(M_D**0.08)*0.453592
    W_ox = ((7*N_pax)**0.702)*0.45359237
    W_fur = 0.412*((N_pax)**1.145)*(((mtow*2.204))**0.489)*0.453592

    weights = np.array([W_wing,W_h,W_v,W_fuselage,W_gear,W_pwr,W_Fuel_sys,W_flightcontrol,W_iae,W_els,W_api,W_ox,W_fur])
    totalweight = sum(weights)

    

    return weights,totalweight
    
def class1est(inputlist,mtowest):
    Cd0 = inputlist[0]
    e = inputlist[1]
    A = inputlist[2]
    propefficiency = inputlist[3]
    u = inputlist[4]
    Vcruise = inputlist[5]
    height = inputlist[6]
    payload = inputlist[7]
    flightrange = inputlist[8]
    loitertime = inputlist[9]

    aref = inputlist[10]
    bref = inputlist[11]

    CLmaxLand = inputlist[12]
    Landingdistance = inputlist[13]

    Cp = u/5918352.

    Vland = np.sqrt(Landingdistance/0.5915)

    WSland = (0.5*CLmaxLand*Vland*Vland*1.225)

    CLcruise = np.sqrt(Cd0*math.pi*e*A)
    CDcruise = 2*Cd0

    CLloiter = np.sqrt(3*Cd0*math.pi*e*A)
    CDloiter = 4*Cd0
    
    LDcruise = CLcruise/CDcruise
    LDloiter = CLloiter/CDloiter

    FfracStart = .99
    FfracTaxi = .995
    FfracTakeoff = .995
    FfracClimb = .985
    FfracDescent = .985

    FfracCruise = 1/math.exp(9.81*Cp*flightrange/(LDcruise*propefficiency))

    WSloiter = WSland*FfracStart*FfracTaxi*FfracTakeoff*FfracClimb*FfracDescent*FfracCruise

    Vloiter = np.sqrt(WSloiter/(CLloiter*0.5*1.225))

    Ffracloiter = 1/math.exp(Vloiter*Cp*loitertime*9.81/(propefficiency*LDloiter))

    FfracTotal = FfracStart*FfracTaxi*FfracTakeoff*FfracClimb*FfracDescent*FfracCruise*Ffracloiter

    Mff =1-FfracTotal

    Mtow = (payload+bref)/(1-aref-Mff)

    oew = aref*Mtow+bref

    fuelw = Mff*Mtow

    return Mtow,oew,fuelw,WSland

def Class2integration(oew,inputlist):

    Cd0 = inputlist[0]
    e = inputlist[1]
    A = inputlist[2]
    propefficiency = inputlist[3]
    u = inputlist[4]
    Vcruise = inputlist[5]
    height = inputlist[6]
    payload = inputlist[7]
    flightrange = inputlist[8]
    loitertime = inputlist[9]

    aref = inputlist[10]
    bref = inputlist[11]

    CLmaxLand = inputlist[12]
    Landingdistance = inputlist[13]

    Cp = u/5918352.

    Vland = np.sqrt(Landingdistance/0.5915)

    WSland = (0.5*CLmaxLand*Vland*Vland*1.225)

    CLcruise = np.sqrt(Cd0*math.pi*e*A)
    CDcruise = 2*Cd0

    CLloiter = np.sqrt(3*Cd0*math.pi*e*A)
    CDloiter = 4*Cd0
    
    LDcruise = CLcruise/CDcruise
    LDloiter = CLloiter/CDloiter

    FfracStart = .99
    FfracTaxi = .995
    FfracTakeoff = .995
    FfracClimb = .985
    FfracDescent = .985

    FfracCruise = 1/math.exp(9.81*Cp*flightrange/(LDcruise*propefficiency))

    WSloiter = WSland*FfracStart*FfracTaxi*FfracTakeoff*FfracClimb*FfracDescent*FfracCruise

    Vloiter = np.sqrt(WSloiter/(CLloiter*0.5*1.225))

    Ffracloiter = 1/math.exp(Vloiter*Cp*loitertime*9.81/(propefficiency*LDloiter))

    FfracTotal = FfracStart*FfracTaxi*FfracTakeoff*FfracClimb*FfracDescent*FfracCruise*Ffracloiter

    Mff =1-FfracTotal

    Mtow = (oew+payload)*1/(1-Mff)

    fuelw = Mff*Mtow

    return Mtow,oew,fuelw

def weightcalcs(inputlist,inputlist2,mtow,oew,fuelw,maxfuelw):
    
    weights,totalweight = class2est(inputlist,inputlist2,mtow,oew,fuelw,maxfuelw)
    mtow,oew,fuelw = Class2integration(totalweight,inputlist)

    return mtow,oew,fuelw

def estweights(inputlist1,inputlist2,deltamaxfuel):

    mtow,oew,fuelw,WSland = class1est(inputlist1,mtowinitialestimate)

    mtowold = mtow
    mtow,oew,fuelw = weightcalcs(inputlist1,inputlist2,mtow,oew,fuelw,fuelw+deltamaxfuel)
    while abs((mtow-mtowold)/mtow) > 0.0001:
        mtowold = mtow
        mtow,oew,fuelw = weightcalcs(inputlist1,inputlist2,mtow,oew,fuelw,fuelw+deltamaxfuel)

    weights,oew = class2est(inputlist,inputlist2,mtow,oew,fuelw,fuelw+deltamaxfuel)

    return mtow,oew,fuelw,WSland,weights
    

Cd0 = 0.02
e = .85
A = 6.1
propefficiency = .82
u = .62
Vcruise = 92.6
height = 18000*0.3048
payload = 444.0
flightrange = 1400000
loitertime = 2700

aref = .614
bref = -5.887

CLmaxLand =2.1
LandingDistance = 500.

mtowinitialestimate = 1700.

sweep = 0.0
sweepquarter = math.pi/180.*sweep
taper = .3
tc = 0.15
Nult = 9.0
rootcoord = 1.76
Maxspeed = 150.0
N_pax = 4
M_D = .4
Sh = 2.24
bh = 2.97
Ah = np.sqrt(bh**2/Sh)
Ch_max = 0.833
tch = .15
trh = tch*Ch_max
lh = 6.0
Sv = 2.0
Av = 3.0
cv_max = 1.34
tcv = 0.15
Trv = tcv*cv_max
bv = np.sqrt(Av*Sv)
fuselage_l = 7.0
fuse_width = 1.2
fuse_height = 1.4
Vc = 92.6
maxfuelw = 400

sep_tanks = 2
We = 130.
Ne = 1

Nult1 = 5.7
lsm = .75
W_L = mtowinitialestimate
integralfraction = 1

inputlist = np.array([Cd0,e,A,propefficiency,u,Vcruise,height,payload,flightrange,loitertime,aref,bref,CLmaxLand,LandingDistance])
inputlist2 = np.array([sweepquarter,taper,tc,Nult,rootcoord,Maxspeed,N_pax,M_D,Sh,Ah,Ch_max,trh,bh,
                       lh,Sv,Av,Trv,bv,cv_max,fuselage_l,fuse_width,fuse_height,Vc,sep_tanks,We,Ne,
                       Nult1,lsm,W_L,integralfraction])

mtow,oew,fuelw,WSland,weights = estweights(inputlist,inputlist2,100)

print "mtow: ", mtow
print "oew:  ", oew
print "fuelw:", fuelw
