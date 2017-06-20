import numpy as np
import VLMtry as VLM
import math
import timeit

###############################################################################
# Aerodynamic center
###############################################################################
# Inputs needed: refx,refy,refz
# need to be defined: wingspan, coord, coords, coordpanels, spanpanels, Force, rho, Vinfx, S 
#raw inputs: refx,refy,refz

def CalcReflength(refx,refy,refz,vortexpoints,wingspan,spanpanels):
    cx = (vortexpoints[:,0,0]-refx)
    cy = (vortexpoints[:,0,1]-refy)+((wingspan/spanpanels)/2)
    cz = (vortexpoints[:,0,2]-refz)
    
    return [cx,cy,cz]
    
def CalcForceCoeff(Force,rho,Vinfx,S):
    # Cx, X-force coefficient ~ Drag
    Cx = Force[:,0]/(rho*Vinfx*Vinfx*S)
    # Cy, Y-force coefficient ~ Side force
    Cy = Force[:,1]/(rho*Vinfx*Vinfx*S)
    # Cz, Z-force coefficient ~ Lift
    Cz = Force[:,2]/(rho*Vinfx*Vinfx*S)

    return [Cx,Cy,Cz]
    
def CalcMomentCoeff(Force,cx,cy,cz,rho,Vinfx,S,wingspan,coord):
    # Cl, rolling moment coefficient
    Mxy = -sum(Force[:,1]*cz)
    Mxz = sum(Force[:,2]*cy)
    l = Mxy+Mxz
    Cl = l/(0.5*rho*Vinfx*Vinfx*S*wingspan)
    # Cm, pitching moment coefficient
    Myx = sum(Force[:,0]*cz)
    Myz = sum(Force[:,2]*cx)
    m = (Myx+Myz)
    Cm = m/(0.5*rho*Vinfx*Vinfx*S*coord)    
    # Cn, Yawing moment coefficient
    Mzx = -sum(Force[:,0]*cy)
    Mzy = sum(Force[:,1]*cx)
    n = Mzx+Mzy
    Cn = n/(0.5*rho*Vinfx*Vinfx*S*wingspan)
    
    return [Cl,Cm,Cn,l,m,n]

def Calcmomentderiv(Force,cx,cy,cz,rho,Vinfx,S,coord,Cl,Cm,Cn,wingspan):
    # dclcy, dCl/dCy
    Mxyly = -sum(((1.01*Force[:,1])*cz))
    Mxzly = sum(Force[:,2]*cy)
    ly = Mxyly+Mxzly
    Cly = ly/(0.5*rho*Vinfx*Vinfx*S*wingspan)
    dclcy = (((Cly-Cl)/Cl)*100)
    # dclcz, dCl/dcz
    Mxylz = -sum(Force[:,1]*cz)
    Mxzlz = sum(((1.01*Force[:,2])*cy))
    lz = Mxylz+Mxzlz
    Clz = lz/(0.5*rho*Vinfx*Vinfx*S*wingspan)
    dclcz = (((Clz-Cl)/Cl)*100)
    # dcmcx, dCm/dCx
    Myxmx = sum(((1.01*Force[:,0])*cz))
    Myzmx = -sum(Force[:,2]*cx)
    mx = Myxmx + Myzmx
    Cmx = mx/(0.5*rho*Vinfx*Vinfx*S*coord)
    dcmcx = (((Cmx-Cm)/Cm)*100)
    # dcmcz, dCm/dcz
    Myxmz = sum(Force[:,0]*cz)
    Myzmz = sum(((1.01*Force[:,2])*cx))
    mz = (Myxmz + Myzmz)
    Cmz = mz/(0.5*rho*Vinfx*Vinfx*S*coord)
    dcmcz = (((Cmz-Cm)/Cm)*100)
    # dcncx, dCn/dCx
    Mzxnx = sum(((1.01*Force[:,0])*cy))
    Mzynx = sum(Force[:,1]*cx)
    nx = Mzxnx+Mzynx
    Cnx = nx/(0.5*rho*Vinfx*Vinfx*S*wingspan)
    dcncx = (((Cnx-Cn)/Cn)*100)
    # dcncy, dCn/dcy
    Mzxny = -sum(Force[:,0]*cy)
    Mzyny = sum(((1.01*Force[:,1])*cx))
    ny = Mzxny+Mzyny
    Cny = ny/(0.5*rho*Vinfx*Vinfx*S*wingspan)
    dcncy = (((Cny-Cn)/Cn)*100)
    
    return [dclcy,dclcz,dcmcx,dcmcz,dcncx,dcncy]

def calcpanelloc(coords,coordpanels,spanpanels,wingspan):
    # calculate the locations of the frontline of the coordpanels
    coordlocx = coords[:-(coordpanels+1),0]
    coordlocx = np.delete(coordlocx,np.arange(coordpanels,coordlocx.size,coordpanels+1))
    # calculate the locations of the (spanwise) midpoint of the spanpanels
    spanlocy = coords[:-(coordpanels+1),1]+(wingspan/(2*spanpanels))
    spanlocy = np.delete(spanlocy,np.arange(coordpanels,spanlocy.size,coordpanels+1))
    # calculate the locations of the height of the panels
    panellocz = coords[:-(coordpanels+1),2]
    panellocz = np.delete(panellocz,np.arange(coordpanels,panellocz.size,coordpanels+1))
    
    return [coordlocx,spanlocy,panellocz]
    
def actualAC(refx,refy,refz,Force,cx,coord):
    # Xac, X-location Aerodynamic center
    FZ1 = Force[:180,2]    
    FZt = sum(FZ1)
    cxred = cx[:180]
    mred = Force[:180,2]*cxred   
    Xac = (sum(mred/FZt))-np.amin(cxred)
    Xacc = Xac/coord
    # Yac, Y-location Aerodynamic center
    Yac = 0
    # Zac, Z-location Aerodyanmic center
    Zac = 0
        
    return [Xac,Yac,Zac,Xacc]

def calcAC(plane,refx,refy,refz,Force,rho,Vinfx):
    
    #General inputs:

    S = plane[10][0][4]
    vortexpoints = plane[4]
    wingspan = plane[10][0][4]
    spanpanels = plane[10][0][6]*6
        
    reflength = CalcReflength(refx,refy,refz,vortexpoints,wingspan,spanpanels)
    
    cx = reflength[0]
    cy = reflength[1]
    cz = reflength[2]

    coord = S/wingspan

    coords = plane[7]
    coordpanels = plane[10][0][5]
    
    forcecoeff = CalcForceCoeff(Force,rho,Vinfx,S)
    momentcoeff = CalcMomentCoeff(Force,cx,cy,cz,rho,Vinfx,S,wingspan,coord)
    
    Cl = momentcoeff[0]
    Cm = momentcoeff[1]
    Cn = momentcoeff[2]
    l = momentcoeff[3]
    m = momentcoeff[4]
    n = momentcoeff[5]
    
    momentderiv = Calcmomentderiv(Force,cx,cy,cz,rho,Vinfx,S,coord,Cl,Cm,Cn,wingspan)
    
    dclcy = momentderiv[0]
    dclcz = momentderiv[1]
    dcmcx = momentderiv[2]
    dcmcz = momentderiv[3]
    dcncx = momentderiv[4]
    dcncy = momentderiv[5]
    
    panellocs = calcpanelloc(coords,coordpanels,spanpanels,wingspan)
    
    coordlocx = panellocs[0]
    spanlocy = panellocs[1]
    panellocz = panellocs[2]
    
    AC = actualAC(refx,refy,refz,Force,m,cx,coord)
    
    return [AC,l,m,n],reflength
    
def calcstabderv1(plane,Vvec,Force,rho,S,reflenght,m,coord):

    cx,cy,cz = reflenght[0],reflenght[1],reflenght[2]
    
    V = np.linalg.norm(Vvec[0])
    Vnew = V+1
    Ratio = Vnew/V

    Vvecv1 = Ratio*Vvec
    Forcev1 = VLM.Calcpoint(plane,Vvecv1,rho)
    Vinfxv1 = Vvecv1[0][0]
    
    dFxv1 = Forcev1[:,0]-Force[:,0]
    dFyv1 = Forcev1[:,1]-Force[:,1]
    dFzv1 = Forcev1[:,2]-Force[:,2]

    dFxu = (Forcev1[:,0]-Force[:,0])
    Cxu = sum(dFxu)/(0.5*rho*Vinfxv1*S)

    dFzu = (Forcev1[:,2]-Force[:,2])
    Czu = sum(dFzu)/(0.5*rho*Vinfxv1*S)

    Myxv1 = sum(Forcev1[:,0]*cz)
    Myzv1 = -sum(Forcev1[:,2]*cx)
    mv1 = Myxv1+Myzv1
    Cmu = (mv1-m)/(0.5*rho*Vinfxv1*Vinfxv1*S*coord)
    
    return Cxu,Czu,Cmu

def calcstabdera1(plane,Vvec,Force,rho,S,reflenght,m,coord):
    
    # Run VLM for a+1 here, outputs should be called Forcea1, Vinfxa1, FLa1, FDa1

    cx,cy,cz = reflenght[0],reflenght[1],reflenght[2]

    alpha = math.atan2(Vvec[0][2],Vvec[0][0])

    deltaalpha = 1

    alpharad = -deltaalpha*math.pi/180.0

    rotmaty = np.matrix([[math.cos(alpharad),0,math.sin(alpharad)],
                     [0,1,0],
                     [-math.sin(alpharad),0,math.cos(alpharad)]])

    alpharad1 = alpha

    alpharad2 = alpha + deltaalpha/180.0*math.pi

    rotmatya1 = np.matrix([[math.cos(alpharad1),0,math.sin(alpharad1)],
                     [0,1,0],
                     [-math.sin(alpharad1),0,math.cos(alpharad1)]])
    rotmatya2 = np.matrix([[math.cos(alpharad2),0,math.sin(alpharad2)],
                     [0,1,0],
                     [-math.sin(alpharad2),0,math.cos(alpharad2)]])

    Vveca1 =np.empty((len(Vvec),3))
    Force1 =np.empty((len(Vvec),3))
    Force2 =np.empty((len(Vvec),3))

    for i in range(len(Vvec)):
        Vveca1[i] = np.dot(rotmaty,Vvec[i])

    Vinfxa1 = Vveca1[0][0]

    Forcea1 = VLM.Calcpoint(plane,Vveca1,rho)

    for i in range(len(Vvec)):
        Force1[i] = np.dot(rotmatya1,Force[i])
        Force2[i] = np.dot(rotmatya2,Forcea1[i])

    FL = sum(Force1[:,2])
    FD = sum(Force1[:,0])

    FLa1 = sum(Force2[:,2])
    FDa1 = sum(Force2[:,0])
    
    dFxa1 = Forcea1[:,0]-Force[:,0]
    dFya1 = Forcea1[:,1]-Force[:,1]
    dFza1 = Forcea1[:,2]-Force[:,2]
    
    CLa = ((FLa1-FL)/(0.5*rho*Vinfxa1*Vinfxa1*S))*57.3
    CDa = ((FDa1-FD)/(0.5*rho*Vinfxa1*Vinfxa1*S))*57.3
    
    Cxa = ((sum(Forcea1[:,0])-sum(Force[:,0]))/(0.5*rho*Vinfxa1*Vinfxa1*S))*57.3
    Cza = ((sum(Forcea1[:,2])-sum(Force[:,2]))/(0.5*rho*Vinfxa1*Vinfxa1*S))*57.3
    
    Myxa1 = sum(Forcea1[:,0]*cz)
    Myza1 = -sum(Forcea1[:,2]*cx)
    ma1 = Myxa1+Myza1
    Cma = ((ma1-m)/(0.5*rho*Vinfxa1*Vinfxa1*S*coord))*57.3

    return CLa,CDa,Cxa,Cza,Cma
    
def calcstabderb1(plane,Vvec,Force,rho,S,reflength,n,l,coord):

    cx,cy,cz = reflength[0],reflength[1],reflength[2]
    
    # Run VLM for b+1 here, outputs should be called Forceb1, Vinfxb1

    Vvecb1 = np.empty((len(Vvec),3))

    betarad = 1.0/180.*math.pi

    rotmatz = np.matrix([[math.cos(betarad),-math.sin(betarad),0],
                         [math.sin(betarad),math.cos(betarad),0],
                         [0,0,1]])
    
    for i in range(len(Vvec)):
        Vvecb1[i] = np.dot(rotmatz,Vvec[i])

    Forceb1 = VLM.Calcpoint(plane,Vvecb1,rho)

    Vinfxb1 = Vvecb1[0][0]
    
    dFxb1 = Forceb1[:,0]-Force[:,0]
    dFyb1 = Forceb1[:,1]-Force[:,1]
    dFzb1 = Forceb1[:,2]-Force[:,2]
    
    Cyb = ((sum(Forceb1[:,1])-sum(Force[:,1]))/(0.5*rho*Vinfxb1*Vinfxb1*S))*57.3
    
    Mxyb1 = -sum(Forceb1[:,1]*cz)
    Myzb1 = sum(Forceb1[:,2]*cy)
    lb1 = Mxyb1+Myzb1
    Clb = ((lb1-l)/(0.5*rho*Vinfxb1*Vinfxb1*S*coord))*57.3
    
    Mzxb1 = sum(Forceb1[:,0]*cy)
    Myzb1 = sum(Forceb1[:,1]*cx)
    nb1 = Mzxb1+Myzb1
    Cnb = ((nb1-n)/(0.5*rho*Vinfxb1*Vinfxb1*S*coord))*57.3
    
    return Cyb,Clb,Cnb

def calcstabderad(deda,CNha,VhV,lhc,ShS):# need CNha, VhV,deda,ShS and lhc either as inputs or pre-defined
    
    Czad = CNha*VhV*VhV*deda*ShS*lhc

    Cmad = -CNha*VhV*VhV*deda*ShS*lhc*lhc
    
    return Czad,Cmad
    
def calcstabderrr(plane,Vvec,Force,rho,S,reflength,n,l,coord,wingspan):

    cx,cy,cz = reflength[0],reflength[1],reflength[2]

    totalpanels = len(plane[5])

    rr = 1.0
    controlpoints = plane[3]
    motionvector = np.array([rr,0,0])

    difvector = np.empty((len(controlpoints),3))

    for i in range(totalpanels):
        difvector[i] = np.cross(controlpoints[i],motionvector)

    # Run VLM for rr+1 here, outputs should be called Forcerr, Vinfxrr

    Vvecrr = Vvec+difvector

    Forcerr = VLM.Calcpoint(plane,Vvecrr,rho)

    Vinfxrr = Vvecrr[0][0]
    
    Mzxrr = -sum(Forcerr[:,0]*cy)
    Mzyrr = sum(Forcerr[:,1]*cx)
    nrr = Mzxrr+Mzyrr
    dnrr = nrr-n
    
    Mxyrr = -sum(Forcerr[:,1]*cz)
    Mxzrr = sum(Forcerr[:,2]*cy)
    lrr = Mxyrr+Mxzrr
    dlrr = lrr-l
    
    dFxrr = Forcerr[:,0]-Force[:,0]
    dFyrr = Forcerr[:,1]-Force[:,1]
    dFzrr = Forcerr[:,2]-Force[:,2]
    
    Cyp = sum(dFyrr)/(0.5*rho*Vinfxrr*Vinfxrr*S)
    
    Clp = (dlrr/(((rr*wingspan)/(2*Vinfxrr))*(0.5*rho*Vinfxrr*Vinfxrr*S*wingspan)))
    Cnp = (dnrr/(((rr*wingspan)/(2*Vinfxrr))*(0.5*rho*Vinfxrr*Vinfxrr*S*wingspan)))
    
    return Cyp,Clp,Cnp

def calcstabderpr(plane,Vvec,Force,rho,S,reflength,m,coord):
    
    # Run VLM for pr+1 here, outputs should be called Forcepr, Vinfxpr

    cx,cy,cz = reflength[0],reflength[1],reflength[2]

    totalpanels = len(plane[5])

    pr = 1.0
    controlpoints = plane[3]
    motionvector = np.array([0,pr,0])

    difvector = np.empty((len(controlpoints),3))

    for i in range(totalpanels):
        difvector[i] = np.cross(controlpoints[i],motionvector)

    Vvecpr1 = Vvec+difvector

    Forcepr = VLM.Calcpoint(plane,Vvecpr1,rho)

    Vinfxpr = Vvec[0][0]

    dFzq = (Forcepr[:,2]-Force[:,2])
    Czq = sum(dFzq)/(((pr*coord)/Vinfxpr)*(0.5*rho*Vinfxpr*Vinfxpr*S*coord)) 
    
    dFxpr = Forcepr[:,0]-Force[:,0]
    dFypr = Forcepr[:,1]-Force[:,1]
    dFzpr = Forcepr[:,2]-Force[:,2]
    
    Myxpr = sum(Forcepr[:,0]*cz)
    Myzpr = -sum(Forcepr[:,2]*cx)
    mpr = Myxpr+Myzpr
    dmpr = mpr-m
    Cmq = (dmpr/(((pr*coord)/Vinfxpr)*(0.5*rho*Vinfxpr*Vinfxpr*S*coord*coord)))
    
    return Czq,Cmq
    
def calcstabderyr(plane,Vvec,Force,rho,S,reflength,n,l,coord,wingspan):            
    # Run VLM for yr+1 here, outputs should be called Forceyr, Vinfxyr

    cx,cy,cz = reflength[0],reflength[1],reflength[2]

    totalpanels = len(plane[5])

    yr = 1.0
    controlpoints = plane[3]
    motionvector = np.array([0,0,yr])

    difvector = np.empty((len(controlpoints),3))

    for i in range(totalpanels):
        difvector[i] = np.cross(controlpoints[i],motionvector)

    Vvecyr1 = Vvec+difvector

    Forceyr = VLM.Calcpoint(plane,Vvecyr1,rho)

    Vinfxyr = Vvec[0][0]

    
    dFxyr = Forceyr[:,0]-Force[:,0]
    dFyyr = Forceyr[:,1]-Force[:,1]
    dFzyr = Forceyr[:,2]-Force[:,2]
    
    Mzxyr = -sum(Forceyr[:,0]*cy)
    Mzyyr = sum(Forceyr[:,1]*cx)
    nyr = Mzxyr+Mzyyr
    dnyr = nyr-n
    
    Mxyyr = -sum(Forceyr[:,1]*cz)
    Mxzyr = sum(Forceyr[:,2]*cy)
    lyr = Mxyyr+Mxzyr
    dlyr = lyr-l
    
    Cl = lyr/(0.5*rho*Vinfxyr*Vinfxyr*S*wingspan)
    
    Clr = (dlyr/(((yr*wingspan)/(2*Vinfxyr))*(0.5*rho*Vinfxyr*Vinfxyr*S*wingspan)))
    Cnr = (dnyr/(((yr*wingspan)/(2*Vinfxyr))*(0.5*rho*Vinfxyr*Vinfxyr*S*wingspan)))
    
    return Clr,Cnr

def stabderives(plane,Vvec,stabinput,rho,deda,CNah,VhV):

    Force = VLM.Calcpoint(plane,Vvec,rho)

    refx = stabinput[0]
    refy = stabinput[1]
    refz = stabinput[2]

    Vinfx = Vvec[0][0]

    AC,reflength = calcAC(plane,refx,refy,refz,Force,rho,Vinfx)
    
    l = AC[1]
    m = AC[2]
    n = AC[3]

    S = plane[10][0][3]
    b = plane[10][0][4]

    lhc = plane[10][1][5]-plane[10][0][7]
    ShS = plane[10][1][2]/S

    coord = S/b
    
    Cxu,Czu,Cmu = calcstabderv1(plane,Vvec,Force,rho,S,reflength,m,coord)
    
    CLa,CDa,Cxa,Cza,Cma = calcstabdera1(plane,Vvec,Force,rho,S,reflength,m,coord)

    Cyb,Clb,Cnb = calcstabderb1(plane,Vvec,Force,rho,S,reflength,n,l,coord)

    Czad,Cmad = calcstabderad(deda,CNah,VhV,lhc,ShS)
    
    Cyp,Clp,Cnp = calcstabderrr(plane,Vvec,Force,rho,S,reflength,n,l,coord,b)
    
    Czq,Cmq = calcstabderpr(plane,Vvec,Force,rho,S,reflength,m,coord)
    
    Clr,Cnr  = calcstabderyr(plane,Vvec,Force,rho,S,reflength,n,l,coord,b)
    
    stabderivs = [Cxu,Czu,Cmu,CLa,CDa,Cxa,Cza,Cma,Cyb,Clb,Cnb,Czad,Cmad,Cyp,Clp,Cnp,Czq,Cmq,Clr,Cnr]    
    
    return stabderivs
    
    
    
    
    
    
    
    
    
    
    
    
    
    

