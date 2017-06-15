import math
import numpy as np
import scipy
import scipy.linalg as la
import timeit

import VLMviscous as VLMV

# utility function to quicky print the coordinates of a given panel in the mesh
def printpanel(a,panellist,coordlist):
    print coordlist[panellist[a][0]]
    print coordlist[panellist[a][1]]
    print coordlist[panellist[a][2]]
    print coordlist[panellist[a][3]]

# heritage function. used to calculate the influece coefficient in the matrix for
# a certain panel on a certain controlpoint
def calcinfluence(panel,controlpoint,vortexpoints,controlpoints,normal):

    # make local copies of the specific vortexpoint
    vtpoint1 = np.array(vortexpoints[panel][0])
    vtpoint2 = np.array(vortexpoints[panel][1])

    # define downstream points for the vortices
    vtpoint1inf = np.array([10000000.0,vtpoint1[1],vtpoint1[2]])
    vtpoint2inf = np.array([10000000.0,vtpoint2[1],vtpoint2[2]])

    # local copy of the control points and normal vector
    ctrpoint = np.array(controlpoints[controlpoint])
    nrml = np.array(normal)
    
    # distance vectors between vortexpoints and controlpoint
    r0 = vtpoint2-vtpoint1
    r1 = ctrpoint-vtpoint1
    r2 = ctrpoint-vtpoint2

    # calculate important parameters for the front vortex, see sources: http://www.dept.aoe.vt.edu/~mason/Mason_f/CAtxtChap6.pdf
    omega1 = np.dot(r0,(r1/np.linalg.norm(r1)))- np.dot(r0,(r2/np.linalg.norm(r2)))
    crossproduct = np.cross(r1,r2)
    psi1 = crossproduct/(np.linalg.norm(crossproduct)**2)
    Vab = psi1*omega1/(4*math.pi)

    # repeat above process for infinite vortexfilament to the rear left
    r0 = vtpoint1-vtpoint1inf
    r1 = ctrpoint-vtpoint1inf
    r2 = ctrpoint-vtpoint1

    omega1 = np.dot(r0,(r1/np.linalg.norm(r1)))- np.dot(r0,(r2/np.linalg.norm(r2)))
    crossproduct = np.cross(r1,r2)
    psi1 = crossproduct/(np.linalg.norm(crossproduct)**2)
    Vainf = psi1*omega1/(4*math.pi)
    
    # repeat process for infinite vortexfilament to the rear right
    r0 = vtpoint2inf-vtpoint2
    r1 = ctrpoint-vtpoint2
    r2 = ctrpoint-vtpoint2inf

    omega1 = np.dot(r0,(r1/np.linalg.norm(r1)))- np.dot(r0,(r2/np.linalg.norm(r2)))
    crossproduct = np.cross(r1,r2)
    psi1 = crossproduct/(np.linalg.norm(crossproduct)**2)
    Vbinf = psi1*omega1/(4*math.pi)

    # calculate total vortex influence
    Vinfl = Vab+Vainf+Vbinf

    # calculate influence coefficient
    infl = np.dot(Vinfl,normal)

    return infl

# calculate the induced velocities based on the geometry
def calcVind(Vmat,vortexpoints,normal):

    # define amount of panels used
    totalpanels = len(normal)
    
    # create local copies of the vortexpoints
    x1 = vortexpoints[:,0,0]
    x2 = vortexpoints[:,1,0]

    y1 = vortexpoints[:,0,1]
    y2 = vortexpoints[:,1,1]

    z1 = vortexpoints[:,0,2]
    z2 = vortexpoints[:,1,2]

    # define point where forces latch to the mesh
    # slight offset of the bounding vortex in order to prevent nan answers
    # 1 micron does not effect results
    x = (x1+x2)/2+0.000001
    y = (y1+y2)/2
    z = (z1+z2)/2

    x2x1 = x2-x1
    y2y1 = y2-y1
    z2z1 = z2-z1
    
    # calculate for all control points the influence coefficients
    # using some numpy magic to speed it up
    # for details see above function and the reference: http://www.dept.aoe.vt.edu/~mason/Mason_f/CAtxtChap6.pdf
    for i in range(len(vortexpoints)):
        
        x1val = vortexpoints[i][0][0]
        x2val = vortexpoints[i][1][0]

        y1val = vortexpoints[i][0][1]
        y2val = vortexpoints[i][1][1]

        z1val = vortexpoints[i][0][2]
        z2val = vortexpoints[i][1][2]

        x1 = np.full(totalpanels,x1val)
        x2 = np.full(totalpanels,x2val)

        y1 = np.full(totalpanels,y1val)
        y2 = np.full(totalpanels,y2val)

        z1 = np.full(totalpanels,z1val)
        z2 = np.full(totalpanels,z2val)

        xx1 = x-x1
        xx2 = x-x2

        yy1 = y-y1
        yy2 = y-y2

        zz1 = z-z1
        zz2 = z-z2

        psi1 = (yy1*zz2-yy2*zz1)
        psi2 = (xx1*zz2-xx2*zz1)
        psi3 = (xx1*yy2-xx2*yy1)

        under2 = psi1*psi1 + psi2*psi2 + psi3*psi3

        psi = np.empty((3,totalpanels))
        psi[0,:] = psi1/under2
        psi[1,:] = psi2/under2
        psi[2,:] = psi3/under2

        length1 = np.sqrt(xx1*xx1+yy1*yy1+zz1*zz1)
        length2 = np.sqrt(xx2*xx2+yy2*yy2+zz2*zz2)

        sigma = (x2x1*xx1+y2y1*yy1+z2z1*zz1)/length1-(x2x1*xx2+y2y1*yy2+z2z1*zz2)/length2

        aftpart1 = 1.0 + xx1/length1
        aftpart2 = 1.0 + xx2/length2

        under = (zz1*zz1+yy1*yy1)
        under2 = (zz2*zz2+yy2*yy2)
        
        Vainfy = zz1/under*aftpart1
        Vainfz = -yy1/under*aftpart1
        
        Vbinfy = -zz2/under2*aftpart2
        Vbinfz = yy2/under2*aftpart2

        Vainf = np.zeros((3,totalpanels))
        Vbinf = np.zeros((3,totalpanels))

        Vainf[1,:] = Vainfy
        Vainf[2,:] = Vainfz
        
        Vbinf[1,:] = Vbinfy
        Vbinf[2,:] = Vbinfz

        sigma[i] = 0

        Vind = (psi*sigma + Vainf + Vbinf)/(4*math.pi)

        Vmat[i,:,0] = Vind[0]
        Vmat[i,:,1] = Vind[1]
        Vmat[i,:,2] = Vind[2]

# calculate influence matrix for calculation purposes
# see influence calculation function for details or the sources: http://www.dept.aoe.vt.edu/~mason/Mason_f/CAtxtChap6.pdf
# also uses numpy magic for speed, no touchy.
def CalcInflMatHorseShoe(A,vortexpoints,controlpoints,normal2):

    normal = normal2.T

    totalpanels = len(controlpoints)
    x1 = vortexpoints[:,0,0]
    x2 = vortexpoints[:,1,0]

    y1 = vortexpoints[:,0,1]
    y2 = vortexpoints[:,1,1]

    z1 = vortexpoints[:,0,2]
    z2 = vortexpoints[:,1,2]
    
    x2x1 = x2-x1
    y2y1 = y2-y1
    z2z1 = z2-z1
    
    x = controlpoints[:,0]
    y = controlpoints[:,1]
    z = controlpoints[:,2]

    for i in range(len(vortexpoints)):

        x1val = vortexpoints[i][0][0]
        x2val = vortexpoints[i][1][0]

        y1val = vortexpoints[i][0][1]
        y2val = vortexpoints[i][1][1]

        z1val = vortexpoints[i][0][2]
        z2val = vortexpoints[i][1][2]

        x1 = np.full(totalpanels,x1val)
        x2 = np.full(totalpanels,x2val)

        y1 = np.full(totalpanels,y1val)
        y2 = np.full(totalpanels,y2val)

        z1 = np.full(totalpanels,z1val)
        z2 = np.full(totalpanels,z2val)

        xx1 = x-x1
        xx2 = x-x2

        yy1 = y-y1
        yy2 = y-y2

        zz1 = z-z1
        zz2 = z-z2

        psi1 = (yy1*zz2-yy2*zz1)
        psi2 = (xx1*zz2-xx2*zz1)
        psi3 = (xx1*yy2-xx2*yy1)

        under2 = psi1*psi1 + psi2*psi2 + psi3*psi3

        psi = np.empty((3,totalpanels))
        psi[0,:] = psi1/under2
        psi[1,:] = psi2/under2
        psi[2,:] = psi3/under2

        length1 = np.sqrt(xx1*xx1+yy1*yy1+zz1*zz1)
        length2 = np.sqrt(xx2*xx2+yy2*yy2+zz2*zz2)

        sigma = (x2x1*xx1+y2y1*yy1+z2z1*zz1)/length1-(x2x1*xx2+y2y1*yy2+z2z1*zz2)/length2

        aftpart1 = 1.0 + xx1/length1
        aftpart2 = 1.0 + xx2/length2

        under = (zz1*zz1+yy1*yy1)
        under2 = (zz2*zz2+yy2*yy2)
        
        Vainfy = zz1/under*aftpart1
        Vainfz = -yy1/under*aftpart1
        
        Vbinfy = -zz2/under2*aftpart2
        Vbinfz = yy2/under2*aftpart2

        Vainf = np.zeros((3,totalpanels))
        Vbinf = np.zeros((3,totalpanels))

        Vainf[1,:] = Vainfy
        Vainf[2,:] = Vainfz
        
        Vbinf[1,:] = Vbinfy
        Vbinf[2,:] = Vbinfz

        Vind = (psi*sigma + Vainf + Vbinf)/(4*math.pi)
        
        A[:,i] = Vind[0,:]*normal[0,:] + Vind[1,:]*normal[1,:] + Vind[2,:]*normal[2,:]

# take the geometry and create the points for the VLM.
# namelely the vortex and control points
# also included are the normal vectors
def Create_Vortex_Points(coords,panels):

    vortexpoints = []
    controlpoints = []
    normal = np.empty((len(panels),3))

    for i in range(len(panels)):

        p1x = coords[panels[i][0]][0]
        p1y = coords[panels[i][0]][1]
        p1z = coords[panels[i][0]][2]

        p2x = coords[panels[i][1]][0]
        p2y = coords[panels[i][1]][1]
        p2z = coords[panels[i][1]][2]

        p3x = coords[panels[i][2]][0]
        p3y = coords[panels[i][2]][1]
        p3z = coords[panels[i][2]][2]

        p4x = coords[panels[i][3]][0]
        p4y = coords[panels[i][3]][1]
        p4z = coords[panels[i][3]][2]
        
        x1 = (p1x+3*p2x)/4
        x2 = (p4x+3*p3x)/4

        y1 = (p1y+3*p2y)/4
        y2 = (p4y+3*p3y)/4

        z1 = (p1z+3*p2z)/4
        z2 = (p4z+3*p3z)/4

        vortexpoints.append([[x1,y1,z1],[x2,y2,z2]])

        x1ctr = (p2x+3*p1x)/4
        x2ctr = (p3x+3*p4x)/4

        y1ctr = (p2y+3*p2y)/4
        y2ctr = (p3y+3*p4y)/4

        z1ctr = (p2z+3*p1z)/4
        z2ctr = (p3z+3*p4z)/4

        x = (x1ctr+x2ctr)/2
        y = (y1ctr+y2ctr)/2 
        z = (z1ctr+z2ctr)/2

        r1 = np.array([p1x-p2x,p1y-p2y,p1z-p2z])
        r2 = np.array([p3x-p2x,p3y-p2y,p3z-p2z])
        r3 = np.array([p3x-p4x,p3y-p4y,p3z-p4z])
        r4 = np.array([p1x-p4x,p1y-p4y,p1z-p4z])

        controlpoints.append([x,y,z])

        norm1 = np.cross(r1,r2)
        norm2 = np.cross(r3,r4)

        normfinal = ((norm1+norm2)/2)

        normal[i] = normfinal/(np.linalg.norm(normfinal))
        
    controlpoints = np.array(controlpoints)
    vortexpoints = np.array(vortexpoints)
    normal = np.array(normal)

    return controlpoints,vortexpoints,normal

# calculate the actual velocities generated on the vortices by the vortices
def Vinduced(Vmat,Gamma,totalpanels):
    
    Vinducedarr = np.zeros((totalpanels,3))

    for i in range(totalpanels):
        Vinducedarr += Vmat[i]*Gamma[i]

    return Vinducedarr

# initilization function. takes the geometry of the mesh generation proccess
# and produces the relevant matrices and prepares those for fast calculation.

def PreparePlane(coords,panels,compstarts,Sfuse):
    totalpanels = len(panels)
    A = np.empty((totalpanels,totalpanels))
    Vindmat = np.zeros((totalpanels,totalpanels,3))
    controlpoints,vortexpoints,normal = Create_Vortex_Points(coords,panels)
    CalcInflMatHorseShoe(A,vortexpoints,controlpoints,normal)
    calcVind(Vindmat,vortexpoints,normal)
    A,piv = la.lu_factor(A)

    return [A,piv,Vindmat,controlpoints,vortexpoints,normal,compstarts,coords,panels,Sfuse]

def Calcpoint(plane,Vvec,rho):
    totalpanels = len(plane[1])
    b = np.empty(totalpanels)
    for i in range(totalpanels):
        b[i] = np.dot(Vvec[i],plane[5][i])
        
    x = la.lu_solve((plane[0],plane[1]),-b)
    
    Vinducedarr = Vinduced(plane[2],x,totalpanels)
    Vinducedarr2 = Vvec + Vinducedarr

    Force = ((np.cross(Vinducedarr2,(plane[4][:,1]-plane[4][:,0]))*rho).T*x).T
    return Force
"""
wingspan = 10.6
coord = 1.4

alpha = 5.
V = 92.6
rho = .7

coordpanels = 10
spanpanels = 100
totalpanels = coordpanels*spanpanels

dy = wingspan/spanpanels
dx = coord/coordpanels
z = 0.0

coords = []
panels = []
Vindmat = np.zeros((totalpanels,totalpanels,3))

for i in range(spanpanels+1):
    for j in range(coordpanels+1):
        coords.append([j*dx,i*dy-wingspan/2,z])

for i in range(spanpanels):
    for j in range(coordpanels):
        panels.append([i*(coordpanels+1)+j+1,
                       i*(coordpanels+1)+j,
                       i*(coordpanels+1)+j+coordpanels+1,
                       i*(coordpanels+1)+j+coordpanels+2])

coords = np.array(coords)
panels = np.array(panels, dtype=int)

Vinfz = V*math.sin(math.pi/180.*alpha)
Vinfy = 0
Vinfx = V*math.cos(math.pi/180.*alpha)

compstarts = [0,0,0,0]

Vvec = np.full((totalpanels,3),np.array([Vinfx,Vinfy,Vinfz]))

plane = PreparePlane(coords,panels,compstarts)
Force = Calcpoint(plane,Vvec)

Force2 = np.empty((totalpanels,3))
alpharad = math.pi/180.*alpha

rotmaty = np.matrix([[math.cos(alpharad),0,math.sin(alpharad)],
                     [0,1,0],
                     [-math.sin(alpharad),0,math.cos(alpharad)]])

time4 = timeit.time.clock()
for i in range(totalpanels):
    Force2[i] = np.dot(rotmaty,Force[i])

Fx = sum(Force[:,0])
Fy = sum(Force[:,1])
Fz = sum(Force[:,2])

FL = sum(Force2[:,2])
FD = sum(Force2[:,0])

S = coord*wingspan
A = wingspan*wingspan/S

CL = FL/(0.5*rho*V*V*S)
CD = FD/(0.5*rho*V*V*S)

e = (CL*CL)/(math.pi*A*CD)
time5 = timeit.time.clock()
"""
