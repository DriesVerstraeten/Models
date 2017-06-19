import VLMmesh as VLMM
import VLMtry as VLM
import VLMviscous as VLMV
import numpy as np
import math
import timeit
import os

import CombinedACStabDerivs as VLMstab

def CreatePlane(wingarray,Hwingarray,Vwingarray,foilpolarpath,foilpaths,tfoil,foilarray,ellipsearray,
                ellipselocarray,PanelsPerHalfHoop,fuselagelength,panelspersection):

    viscdata = VLMV.importpolars(foilpolarpath)

    planedata = [wingarray,Hwingarray,Vwingarray,foilpolarpath,foilpaths,tfoil,foilarray,ellipsearray,
                ellipselocarray,PanelsPerHalfHoop,fuselagelength,panelspersection]

    coords,uppers,lowers = VLMM.loadfoils(foilpaths)
    
    planecoords,planepanels,Sfuse = VLMM.CreatePlaneGeom(wingarray,Hwingarray,Vwingarray,ellipsearray,ellipselocarray,PanelsPerHalfHoop,
                                                         fuselagelength,panelspersection,foilarray,foilpaths,tfoil)

    tailcoords,tailpanels = VLMM.createHtailpoints(Hwingarray[0],Hwingarray[1],Hwingarray[2],Hwingarray[3],Hwingarray[4],
                                                   Hwingarray[5],Hwingarray[6],Hwingarray[7],Hwingarray[8],tfoil,coords,uppers,lowers)
    #tailcoords,tailpanels = VLMM.createHtailpoints(taperh,bh,Sh,xpanelsh,SpanPanelsPerSectionh,xwingh,fuseelipsewidthh,zwingh,dihedralh,tfoil,coords,uppers,lowers)

    plane = VLM.PreparePlane(planecoords,planepanels,[comp1,comp2,comp3,comp4],Sfuse)
    tail = VLM.PreparePlane(tailcoords,tailpanels,[comp1,comp2,comp3,comp4],Sfuse)

    plane.append([wingarray,Hwingarray,Vwingarray,foilpolarpath,foilpaths,tfoil,foilarray,ellipsearray,
                ellipselocarray,PanelsPerHalfHoop,fuselagelength,panelspersection])

    tail.append([wingarray,Hwingarray,Vwingarray,foilpolarpath,foilpaths,tfoil,foilarray,ellipsearray,
                ellipselocarray,PanelsPerHalfHoop,fuselagelength,panelspersection])

    return plane,tail,viscdata

def CalcTailPars(plane,tail,data):

    V = 100.0
    rho = 1.225

    Sref = plane[10][0][3]

    Sreftail = plane[10][1][2]

    alpha = 1.0
    alpha2 = 2.0

    tailstart = plane[6][1]
    tailend = plane[6][2]

    totalpanels = len(plane[5])
    totaltailpanels = len(tail[5])

    Vinfz = V*math.sin(math.pi/180.*alpha)
    Vinfy = 0
    Vinfx = V*math.cos(math.pi/180.*alpha)

    Vinfz2 = V*math.sin(math.pi/180.*alpha2)
    Vinfy2 = 0
    Vinfx2 = V*math.cos(math.pi/180.*alpha2)

    Vvec1 = np.full((totalpanels,3),np.array([Vinfx,Vinfy,Vinfz]))
    Vvec2 = np.full((totalpanels,3),np.array([Vinfx2,Vinfy2,Vinfz2]))

    Vvectail1 = Vvec1[:totaltailpanels]
    Vvectail2 = Vvec2[:totaltailpanels]
    
    Forces1 = VLM.Calcpoint(plane,Vvec1,rho)
    TailForces1 = VLM.Calcpoint(tail,Vvectail1,rho)

    Forces2 = VLM.Calcpoint(plane,Vvec2,rho)
    TailForces2 = VLM.Calcpoint(tail,Vvectail2,rho)

    alpharad = math.pi/180.*alpha

    rotmaty = np.matrix([[math.cos(alpharad),0,math.sin(alpharad)],
                         [0,1,0],
                         [-math.sin(alpharad),0,math.cos(alpharad)]])

    alpharad2 = math.pi/180.*alpha2

    rotmaty2 = np.matrix([[math.cos(alpharad2),0,math.sin(alpharad2)],
                         [0,1,0],
                         [-math.sin(alpharad2),0,math.cos(alpharad2)]])

    AeroFrameF1 = np.empty((totalpanels,3))
    AeroFrameTF1 = np.empty((totaltailpanels,3))

    AeroFrameF2 = np.empty((totalpanels,3))
    AeroFrameTF2 =np.empty((totaltailpanels,3))

    for i in range(totalpanels):
        AeroFrameF1[i] = np.dot(rotmaty,Forces1[i])
        AeroFrameF2[i] = np.dot(rotmaty2,Forces2[i])

    for i in range(totaltailpanels):
        AeroFrameTF1[i] = np.dot(rotmaty,TailForces1[i])
        AeroFrameTF2[i] = np.dot(rotmaty2,TailForces2[i])

    FL1 = sum(AeroFrameF1[:,2])
    FL2 = sum(AeroFrameF2[:,2])

    TailFL1 = sum(AeroFrameF1[tailstart:tailend,2])
    TailFL2 = sum(AeroFrameF2[tailstart:tailend,2])

    TailFL1Clean = sum(AeroFrameTF1[:,2])
    TailFL2Clean = sum(AeroFrameTF2[:,2])

    CL1 = FL1/(0.5*rho*V*V*Sref)
    CL2 = FL2/(0.5*rho*V*V*Sref)

    CLtail1 = TailFL1/(0.5*rho*V*V*Sreftail)
    CLtail2 = TailFL2/(0.5*rho*V*V*Sreftail)

    CLminH1 = (FL1-TailFL1)/(0.5*rho*V*V*Sref)
    CLminH2 = (FL2-TailFL2)/(0.5*rho*V*V*Sref)

    CLalphaPlane =  (CL2-CL1)*180/math.pi
    CLalphaminH = (CLminH2-CLminH1)*180/math.pi

    CLH1 = TailFL1Clean/(0.5*rho*V*V*Sreftail)
    CLH2 = TailFL2Clean/(0.5*rho*V*V*Sreftail)

    CLalphaH = (CLH2-CLH1)*180/math.pi
    CLalphaHdirty = (CLtail2-CLtail1)*180/math.pi

    deda = CLalphaHdirty/CLalphaH

    return deda,CLalphaminH,CLalphaH

def CalcCLCD(plane,V,alpha,rho,data,DynViscosity):

    totalpanels = len(plane[5])

    alpharad = np.radians(alpha)
    
    Vinfz = V*math.sin(math.pi/180.*alpha)
    Vinfy = 0
    Vinfx = V*math.cos(math.pi/180.*alpha)

    Vvec = np.full((totalpanels,3),np.array([Vinfx,Vinfy,Vinfz]))

    Force = VLM.Calcpoint(plane,Vvec,rho)

    Force2 = np.empty((totalpanels,3))

    rotmaty = np.matrix([[math.cos(alpharad),0,math.sin(alpharad)],
                         [0,1,0],
                         [-math.sin(alpharad),0,math.cos(alpharad)]])

    for i in range(totalpanels):
        Force2[i] = np.dot(rotmaty,Force[i])

    FL = sum(Force2[:,2])
    FD = sum(Force2[:,0])

    wingarray = plane[10][0]
    Hwingarray = plane[10][1]
    Vwingarray = plane[10][2]
    foilarray = plane[10][6]
    tfoil = plane[10][5]
    fuselagelength = plane[10][10]

    VDrag = VLMV.ViscousDrag(plane,Force,wingarray,Hwingarray,Vwingarray,data,V,rho,foilarray,tfoil,DynViscosity,plane[9],alpha,fuselagelength)

    return FL,FD,VDrag

def CreateCLCDpolar(plane,data,V,rho,DynViscosity,alphastart,alphaend,alphapoints):

    points = []
    
    for alpha in np.linspace(alphastart,alphaend,alphapoints):
        
        FL,FD,FDv = CalcCLCD(plane,V,alpha,rho,data,DynViscosity)

        points.append([alpha,FL,FD,FDv,FD+FDv])

    points = np.array(points)

    return points

def CalcFL(plane,alpha,V,rho):

    Vinfz = V*math.sin(math.pi/180.*alpha)
    Vinfy = 0
    Vinfx = V*math.cos(math.pi/180.*alpha)

    Vvec = np.full((totalpanels,3),np.array([Vinfx,Vinfy,Vinfz]))

    Force = VLM.Calcpoint(plane,Vvec,rho)

    Force2 = np.empty((totalpanels,3))
        
    alpharad = math.pi/180.*alpha

    rotmaty = np.matrix([[math.cos(alpharad),0,math.sin(alpharad)],
                         [0,1,0],
                         [-math.sin(alpharad),0,math.cos(alpharad)]])
    
    for i in range(totalpanels):
        Force2[i] = np.dot(rotmaty,Force[i])

    return Force2
    

def CalcDragForConddition(plane,W,V,rho,DynViscosity,data):

    alphanew = 0.5
    i = 0
    FLnew = sum(CalcFL(plane,alphanew,V,rho)[:,2])

    while abs(FLnew-W) > .01:
        
        alpha0 = alphanew-0.001
        alpha1 = alphanew+0.001

        FL0 = sum(CalcFL(plane,alpha0,V,rho)[:,2])
        FL1 = sum(CalcFL(plane,alpha1,V,rho)[:,2])

        der = (alpha1-alpha0)/(FL1-FL0)

        alphanew = der*(W-FL0)+alpha0

        FLnew = sum(CalcFL(plane,alphanew,V,rho)[:,2])
        i += 1


    Force = CalcFL(plane,alphanew,V,rho)

    wingarray = plane[10][0]
    Hwingarray = plane[10][1]
    Vwingarray = plane[10][2]
    foilarray = plane[10][6]
    tfoil = plane[10][5]
    fuselagelength = plane[10][10]

    VDrag = VLMV.ViscousDrag(plane,Force,wingarray,Hwingarray,Vwingarray,data,V,rho,foilarray,tfoil,DynViscosity,plane[9],alpha,fuselagelength)

    FD = sum(Force[:,0])

    totalDrag = VDrag + FD

    return totalDrag
    
currpath = os.path.dirname(os.path.abspath(__file__))

foilarray = [0,0,0,0,0,0,0,0]
foilpaths = [currpath+"/foilpolars/finalfoil mod1.dat",currpath+ "/foilpolars/foil1 modified.dat",
             currpath+"/foilpolars/NACA0012.dat"]
foilpolarpath =[currpath+"/foilpolars/finalfoil mod1",currpath+ "/foilpolars/foil1 modified",
                currpath+"/foilpolars/NACA0012"]

tfoil = 1
r2 = .95
r3 = .85
r4 = .4
S = 15.0
b = 10.6
xpanels = 5
SpanPanelsPerSection = 6
xwing = 1.5
fuseelipsewidth = .55
zwing = -0.2
dihedral = 8.

wingarray = [r2,r3,r4,S,b,xpanels,SpanPanelsPerSection,xwing,fuseelipsewidth,zwing,dihedral]

tfoil = 2

taperh = .3
bh = 3.0
Sh = 2.0
xpanelsh = 3
SpanPanelsPerSectionh = 8
xwingh = 8.0
fuseelipsewidthh = 0.15
zwingh = .8
dihedralh = 5.0

Hwingarray = [taperh,bh,Sh,xpanelsh,SpanPanelsPerSectionh,xwingh,fuseelipsewidthh,zwingh,dihedralh]

taperv = 1
bv = 2.0
Sv = 2.0
xpanelsv = 5
SpanPanelsPerSectionv = 8
xwingv = 8.0
fuseelipsewidthv = .3
zwingv = .6

comp1 = 0
comp2 = xpanels*SpanPanelsPerSection*6
comp3 = comp2 + xpanelsh*SpanPanelsPerSectionh*2
comp4 = comp3 + xpanelsv*SpanPanelsPerSectionv

Vwingarray = [taperv,bv,Sv,xpanelsv,SpanPanelsPerSectionv,xwingv,fuseelipsewidthv,zwingv]

hoop1 = np.array([0.3,0.2,0.2])
hoop1loc = np.array([0.2,0.0,0.3])

hoop2 = np.array([0.5,0.3,0.3])
hoop2loc = np.array([1.0,0.0,0.3])

hoop3 = np.array([0.6,0.4,0.6])
hoop3loc = np.array([2.0,0.0,0.2])

hoop4 = np.array([0.6,0.4,1.0])
hoop4loc = np.array([3.0,0.0,0.2])

hoop5 = np.array([0.6,0.4,.6])
hoop5loc = np.array([4.0,0.0,0.3])

hoop6 = np.array([0.4,0.3,.4])
hoop6loc = np.array([5.0,0.0,0.4])

hoop7 = np.array([0.2,0.2,.2])
hoop7loc = np.array([6.0,0.0,0.55])

ellipsearray = np.array([hoop1,hoop2,hoop3,hoop4,hoop5,hoop6,hoop7])
ellipselocarray = np.array([hoop1loc,hoop2loc,hoop3loc,hoop4loc,hoop5loc,hoop6loc,hoop7loc])

PanelsPerHalfHoop = 6
fuselagelength = 7.0
panelspersection = [1,1,1,1,1,1,1,1]


plane,tail,viscdata = CreatePlane(wingarray,Hwingarray,Vwingarray,foilpolarpath,foilpaths,tfoil,foilarray,ellipsearray,
                ellipselocarray,PanelsPerHalfHoop,fuselagelength,panelspersection)


deda,CNAminH,CNah = CalcTailPars(plane,tail,viscdata)
VhV = 1.0

refx = 0.0
refy = 0.0
refz = 0.0

stabinputs = [refx,refy,refz]

DynViscosity = 1.725*10**-5

alpha = 4.
V = 92.6
rho = 1.225

alphastart = 0.0
alphaend = 5.0
alphapoints = 10

totalpanels = len(plane[5])

Vinfz = V*math.sin(math.pi/180.*alpha)
Vinfy = 0
Vinfx = V*math.cos(math.pi/180.*alpha)

Vvec = np.full((totalpanels,3),np.array([Vinfx,Vinfy,Vinfz]))

stabderives = VLMstab.stabderives(plane,Vvec,stabinputs,rho,deda,CNah,VhV)

