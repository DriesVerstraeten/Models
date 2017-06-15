import VLMmesh as VLMM
import VLMtry as VLM
import VLMviscous as VLMV
import numpy as np
import math

foilarray = [0,0,0,0,0,0,0,0]
foilpaths = ["C:/Users/Jaep-PC/Documents/foilpolars/finalfoil mod1.dat","C:/Users/Jaep-PC/Documents/foilpolars/foil1 modified.dat",
             "C:/Users/Jaep-PC/Documents/foilpolars/NACA0012.dat"]
foilpolarpath =["C:/Users/Jaep-PC/Documents/foilpolars/finalfoil mod1","C:/Users/Jaep-PC/Documents/foilpolars/foil1 modified",
                "C:/Users/Jaep-PC/Documents/foilpolars/NACA0012"]
tfoil = 1
r2 = .95
r3 = .85
r4 = .4
S = 15.0
b = 10.6
xpanels = 5
SpanPanelsPerSection = 6
xwing = 1.5
fuseelipsewidth = 0.02
zwing = -0.2
dihedral = 8.

wingarray = [r2,r3,r4,S,b,xpanels,SpanPanelsPerSection,xwing,fuseelipsewidth,zwing,dihedral]

coords,uppers,lowers = VLMM.loadfoils(foilpaths)

tfoil = 2

taperh = .3
bh = 3.5
Sh = 2.0
xpanelsh = 3
SpanPanelsPerSectionh = 8
xwingh = 8.0
fuseelipsewidthh = 1.1
zwingh = .6
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

planecoords,planepanels,Sfuse = VLMM.CreatePlaneGeom(wingarray,Hwingarray,Vwingarray,ellipsearray,ellipselocarray,PanelsPerHalfHoop,fuselagelength,panelspersection,foilarray,foilpaths,tfoil)
tailcoords,tailpanels = VLMM.createHtailpoints(taperh,bh,Sh,xpanelsh,SpanPanelsPerSectionh,xwingh,fuseelipsewidthh,zwingh,dihedralh,tfoil,coords,uppers,lowers)

totalpanels = len(planepanels)
totaltailpanels = len(tailpanels)

alpha = 3.0
V = 92.6
rho = .7

Vinfz = V*math.sin(math.pi/180.*alpha)
Vinfy = 0.0
Vinfx = V*math.cos(math.pi/180.*alpha)
DynViscosity = 1.725*10**-5

Vvec = np.full((totalpanels,3),np.array([Vinfx,Vinfy,Vinfz]))
Vvectail = np.full((totaltailpanels,3),np.array([Vinfx,Vinfy,Vinfz]))

data = VLMV.importpolars(foilpolarpath)

plane = VLM.PreparePlane(planecoords,planepanels,[comp1,comp2,comp3,comp4],Sfuse)
tail = VLM.PreparePlane(tailcoords,tailpanels,[comp1,comp2,comp3,comp4],Sfuse)

Forces = VLM.Calcpoint(plane,Vvec,rho)
TailForces = VLM.Calcpoint(tail,Vvectail,rho)

VDrag = VLMV.ViscousDrag(plane,Forces,wingarray,Hwingarray,Vwingarray,data,V,rho,foilarray,tfoil,DynViscosity,plane[9],alpha,fuselagelength)
