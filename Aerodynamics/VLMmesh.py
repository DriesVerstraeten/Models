import math
import numpy as np
import matplotlib.pyplot as plt
import timeit


def Interpolate3D(pointa,pointb,c):

    pointc = pointa+(pointb-pointa)*c

    return pointc

def GetSurfaceY(surface,x):

    i = 0

    dif = x-surface[i][0]

    if dif > 0:
        while x > surface[i][0]:
            i+=1
    if dif < 0:
        while x < surface[i][0]:
            i+=1

    i -= 1

    x0 = surface[i][0]
    x1 = surface[i+1][0]

    y0 = surface[i][1]
    y1 = surface[i+1][1]

    y = y0 + (x-x0)*(y1-y0)/(x1-x0)

    return y

def GetCamberY(upper,lower,x):

    return (GetSurfaceY(upper,x)+GetSurfaceY(lower,x))/2

def FoilCoordDef(foil):

    a = np.genfromtxt(foil,skip_header = 1)

    uppersurface = []
    lowersurface = []

    uppersurface.append(a[0])

    switch = False

    for i in range(1,len(a)):
        if(a[i-1][0] > a[i][0]):
            uppersurface.append(a[i])
        else:
            if switch ==False:
                lowersurface.append(a[i-1])
                switch = True
            lowersurface.append(a[i])

    uppersurface = np.array(uppersurface)
    lowersurface = np.array(lowersurface)

    return a,uppersurface,lowersurface

def ChordBase(r2,r3,r4,b,S):

    return S/((r4/2.0+r3+r2+.5)*(b/3.0))

def GetChords(r2,r3,r4,b,S):

    C1 = ChordBase(r2,r3,r4,b,S)
    C2 = C1*r2
    C3 = C1*r3
    C4 = C1*r4

    return C1,C2,C3,C4

def createwingpoints(r2,r3,r4,b,S,xpanels,SpanPanelsPerSection,xwing,fuseelipsewidth,zwing,dihedral,foilarray,coords,uppers,lowers):
    C1,C2,C3,C4 = GetChords(r2,r3,r4,b,S)

    C = np.array([C1,C2,C3,C4])

    StartingPoints = []

    for i in range(3):
        StartingPoints.append(np.array([xwing-0.25*C[-i+3],-b/2.+b/6.0*i,0]))

    Trapezoidpoint = np.array([xwing-0.25*C[0],0,0])

    c = 1-(fuseelipsewidth/2)/(b/6.0)
    FuseStartPoint = Interpolate3D(StartingPoints[2],Trapezoidpoint,c)
    StartingPoints.append(FuseStartPoint)

    FuseStartPoint2 = np.copy(FuseStartPoint)
    FuseStartPoint2[1] = -FuseStartPoint2[1]
    StartingPoints.append(FuseStartPoint2)

    for i in range(3):
        StartingPoints.append(np.array([xwing-0.25*C[i+1],b/6.0*(i+1),0]))

    StartingPoints = np.array(StartingPoints)

    foilpoints = []
    dx = 1./xpanels
    for j in range(8):

        foilpointslocal = []
        
        for i in range(1,xpanels+1):

            x = dx*i
            y = GetCamberY(uppers[foilarray[j]],lowers[foilarray[j]],x)

            point = np.array([x,0,y])

            foilpointslocal.append(point)
        foilpoints.append(foilpointslocal)

    foilpoints = np.array(foilpoints)
    camberwingpoints = []
    dihedralrad = np.radians(dihedral)

    for i in range(8):

        localfoilpoints = []
        
        if i < 4:
            Cuse = C[-i+3]
        else:
            Cuse = C[i-4]
                
        for j in range(xpanels+1):

            if j == 0:
                localfoilpoints.append(StartingPoints[i])
            else:
                localfoilpoints.append(StartingPoints[i]+Cuse*foilpoints[i][j-1])

        camberwingpoints.append(localfoilpoints)
        
    camberwingpoints = np.array(camberwingpoints)

    rotmatxleftwing =np.matrix([[1,0,0],
                                [0,math.cos(-dihedralrad),-math.sin(-dihedralrad)],
                                [0,math.sin(-dihedralrad),math.cos(-dihedralrad)]])

    rotmatxrightwing = np.matrix([[1,0,0],
                                  [0,math.cos(dihedralrad),-math.sin(dihedralrad)],
                                  [0,math.sin(dihedralrad),math.cos(dihedralrad)]])

    camberwingpointsdih = np.empty(camberwingpoints.shape)

    for i in range(4):
        for j in range(xpanels+1):
            camberwingpointsdih[i][j] = np.dot(rotmatxleftwing,camberwingpoints[i][j])

    for i in range(4,8):
        for j in range(xpanels+1):
            camberwingpointsdih[i][j] = np.dot(rotmatxrightwing,camberwingpoints[i][j])

    camberwingpointsdih[:,:,2] += zwing
    
    dc = 1./(SpanPanelsPerSection)
    wingpoints = []
    
    for i in range(3):
        for j in range(SpanPanelsPerSection):

            c = dc*j
            
            for k in range(xpanels+1):
                wingpoints.append(Interpolate3D(camberwingpointsdih[i][k],camberwingpointsdih[i+1][k],c))
    for i in range(xpanels+1):
        wingpoints.append(camberwingpointsdih[3][i])


    for i in range(4,7):
        for j in range(SpanPanelsPerSection):

            c = dc*j
            
            for k in range(xpanels+1):
                wingpoints.append(Interpolate3D(camberwingpointsdih[i][k],camberwingpointsdih[i+1][k],c))
    for i in range(xpanels+1):
        wingpoints.append(camberwingpointsdih[7][i])

    wingpoints = np.array(wingpoints)

    panels = []

    spanpanels = SpanPanelsPerSection*6
    coordpanels = xpanels

    for i in range(spanpanels/2):
        for j in range(coordpanels):
            panels.append([i*(coordpanels+1)+j+1,
                           i*(coordpanels+1)+j,
                           i*(coordpanels+1)+j+coordpanels+1,
                           i*(coordpanels+1)+j+coordpanels+2])

    for i in range(spanpanels/2+1,spanpanels+1):
        for j in range(coordpanels):
            panels.append([i*(coordpanels+1)+j+1,
                           i*(coordpanels+1)+j,
                           i*(coordpanels+1)+j+coordpanels+1,
                           i*(coordpanels+1)+j+coordpanels+2])

    panels = np.array(panels,dtype = int)

    return wingpoints,panels

def createHtailpoints(taper,bh,S,xpanels,SpanPanelsPerSection,xwing,fuseelipsewidth,zwing,dihedral,foil,coords,uppers,lowers):

    StartingPoints = []

    C0 = 2*S/((1+taper)*bh)
    C1 = taper*C0

    C = np.array([C0,C1])

    StartingPoints.append(np.array([xwing-0.25*C[1],-bh/2.0,0]))

    Trapezoidpoint = np.array([xwing-0.25*C[0],0,0])

    c = 1-(fuseelipsewidth/2)/(bh/2.0)
    FuseStartPoint = Interpolate3D(StartingPoints[0],Trapezoidpoint,c)
    StartingPoints.append(FuseStartPoint)

    FuseStartPoint2 = np.copy(FuseStartPoint)
    FuseStartPoint2[1] = -FuseStartPoint2[1]
    StartingPoints.append(FuseStartPoint2)

    StartingPoints.append(np.array([xwing-0.25*C[1],bh/2.,0]))

    StartingPoints = np.array(StartingPoints)

    foilpoints = []
    dx = 1./xpanels
    for j in range(4):

        foilpointslocal = []
        
        for i in range(1,xpanels+1):

            x = dx*i
            y = GetCamberY(uppers[foil],lowers[foil],x)

            point = np.array([x,0,y])

            foilpointslocal.append(point)
        foilpoints.append(foilpointslocal)

    foilpoints = np.array(foilpoints)
    camberwingpoints = []
    dihedralrad = np.radians(dihedral)

    for i in range(4):

        localfoilpoints = []
        
        if i < 2:
            Cuse = C[-i+1]
        else:
            Cuse = C[i-2]
                
        for j in range(xpanels+1):

            if j == 0:
                localfoilpoints.append(StartingPoints[i])
            else:
                localfoilpoints.append(StartingPoints[i]+Cuse*foilpoints[i][j-1])

        camberwingpoints.append(localfoilpoints)
        
    camberwingpoints = np.array(camberwingpoints)

    rotmatxleftwing =np.matrix([[1,0,0],
                                [0,math.cos(-dihedralrad),-math.sin(-dihedralrad)],
                                [0,math.sin(-dihedralrad),math.cos(-dihedralrad)]])

    rotmatxrightwing = np.matrix([[1,0,0],
                                  [0,math.cos(dihedralrad),-math.sin(dihedralrad)],
                                  [0,math.sin(dihedralrad),math.cos(dihedralrad)]])

    camberwingpointsdih = np.empty(camberwingpoints.shape)

    for i in range(2):
        for j in range(xpanels+1):
            camberwingpointsdih[i][j] = np.dot(rotmatxleftwing,camberwingpoints[i][j])

    for i in range(2,4):
        for j in range(xpanels+1):
            camberwingpointsdih[i][j] = np.dot(rotmatxrightwing,camberwingpoints[i][j])

    camberwingpointsdih[:,:,2] += zwing
    
    dc = 1./(SpanPanelsPerSection)
    wingpoints = []
    
    for i in range(1):
        for j in range(SpanPanelsPerSection):

            c = dc*j
            
            for k in range(xpanels+1):
                wingpoints.append(Interpolate3D(camberwingpointsdih[i][k],camberwingpointsdih[i+1][k],c))
    for i in range(xpanels+1):
        wingpoints.append(camberwingpointsdih[1][i])


    for i in range(2,3):
        for j in range(SpanPanelsPerSection):

            c = dc*j
            
            for k in range(xpanels+1):
                wingpoints.append(Interpolate3D(camberwingpointsdih[i][k],camberwingpointsdih[i+1][k],c))
    for i in range(xpanels+1):
        wingpoints.append(camberwingpointsdih[3][i])

    wingpoints = np.array(wingpoints)

    panels = []

    spanpanels = SpanPanelsPerSection*2
    coordpanels = xpanels

    for i in range(spanpanels/2):
        for j in range(coordpanels):
            panels.append([i*(coordpanels+1)+j+1,
                           i*(coordpanels+1)+j,
                           i*(coordpanels+1)+j+coordpanels+1,
                           i*(coordpanels+1)+j+coordpanels+2])

    for i in range(spanpanels/2+1,spanpanels+1):
        for j in range(coordpanels):
            panels.append([i*(coordpanels+1)+j+1,
                           i*(coordpanels+1)+j,
                           i*(coordpanels+1)+j+coordpanels+1,
                           i*(coordpanels+1)+j+coordpanels+2])

    panels = np.array(panels,dtype = int)

    return wingpoints,panels

def createVtailpoints(taper,bv,S,xpanels,SpanPanelsPerSection,xwing,fuseelipsewidth,zwing,foil,coords,uppers,lowers):

    StartingPoints = []

    C0 = 2*S/((1+taper)*bv)
    C1 = taper*C0

    C = np.array([C0,C1])

    StartingPoints.append(np.array([xwing-0.25*C[1],0,bv]))

    Trapezoidpoint = np.array([xwing-0.25*C[0],0,0])

    c = 1-(fuseelipsewidth/2)/(bv)
    FuseStartPoint = Interpolate3D(StartingPoints[0],Trapezoidpoint,c)
    StartingPoints.append(FuseStartPoint)

    StartingPoints = np.array(StartingPoints)

    foilpoints = []
    dx = 1./xpanels
    for j in range(2):

        foilpointslocal = []
        
        for i in range(1,xpanels+1):

            x = dx*i
            y = GetCamberY(uppers[foil],lowers[foil],x)

            point = np.array([x,y,0])

            foilpointslocal.append(point)
        foilpoints.append(foilpointslocal)

    foilpoints = np.array(foilpoints)
    camberwingpoints = []

    for i in range(2):

        localfoilpoints = []
                
        for j in range(xpanels+1):

            if j == 0:
                localfoilpoints.append(StartingPoints[i])
            else:
                localfoilpoints.append(StartingPoints[i]+C[i-1]*foilpoints[i][j-1])

        camberwingpoints.append(localfoilpoints)
        
    camberwingpointsdih = np.array(camberwingpoints)

    camberwingpointsdih[:,:,2] += zwing
    
    dc = 1./(SpanPanelsPerSection)
    wingpoints = []
    
    for i in range(1):
        for j in range(SpanPanelsPerSection):

            c = dc*j
            
            for k in range(xpanels+1):
                wingpoints.append(Interpolate3D(camberwingpointsdih[i][k],camberwingpointsdih[i+1][k],c))
    for i in range(xpanels+1):
        wingpoints.append(camberwingpointsdih[1][i])

    wingpoints = np.array(wingpoints)

    panels = []

    spanpanels = SpanPanelsPerSection*2
    coordpanels = xpanels

    for i in range(spanpanels/2):
        for j in range(coordpanels):
            panels.append([i*(coordpanels+1)+j+1,
                           i*(coordpanels+1)+j,
                           i*(coordpanels+1)+j+coordpanels+1,
                           i*(coordpanels+1)+j+coordpanels+2])

    panels = np.array(panels,dtype = int)

    return wingpoints,panels

def loadfoils(patharray):

    coordarray = []
    upperarray = []
    lowerarray = []

    for i in range(len(patharray)):
        coords,upper,lower = FoilCoordDef(patharray[i])

        coordarray.append(coords)
        upperarray.append(upper)
        lowerarray.append(lower)

    return coordarray,upperarray,lowerarray

def CreateHoopPoints(a,b1,b2,PanelsPerHalfHoop):

    hooppoints = []

    #lower hoops

    for i in np.linspace(0.0,math.pi,PanelsPerHalfHoop+1):

        y = a*math.cos(i)
        z = -np.sqrt(b1*b1*(1-y*y/(a*a)))
        x = 0
        
        hooppoints.append(np.array([x,y,z]))

    #upperhoops
    j = 0
    for i in np.linspace(0.0,math.pi,PanelsPerHalfHoop+1):

        y = -a*math.cos(i)
        z = np.sqrt(b2*b2*(1-y*y/(a*a)))
        x = 0
        if (j != 0) and (j != PanelsPerHalfHoop):
            hooppoints.append(np.array([x,y,z]))

        j+=1

    hooppoints = np.array(hooppoints)

    return hooppoints


def CreateFuseMesh(ellipsearray,ellipselocarray,PanelsPerHalfHoop,PanelsPerSectionlist,fuselagelength):

    hooppointlist = []

    for i in range(len(ellipsearray)):

        hooppos = ellipselocarray[i]

        a = ellipsearray[i][0]
        b1 = ellipsearray[i][1]
        b2 = ellipsearray[i][2]

        hooppointlist.append(hooppos+CreateHoopPoints(a,b1,b2,PanelsPerHalfHoop))
        
    hooppointlist = np.array(hooppointlist)
    frontpoint = np.copy(ellipselocarray[0])
    aftpoint = np.copy(ellipselocarray[-1])
    
    frontpoint[0] = 0
    aftpoint[0] = fuselagelength

    FusePoints = []

    FusePoints.append(frontpoint)
    dc = 1./PanelsPerSectionlist[0]

    for i in range(PanelsPerSectionlist[0]):

        c = (i+1)*dc
        for j in range(PanelsPerHalfHoop*2): 
            FusePoints.append(Interpolate3D(frontpoint,hooppointlist[0][j],c))
    
    for i in range(0,len(ellipselocarray)-1):
        dc = 1./PanelsPerSectionlist[i+1]

        for j in range(PanelsPerSectionlist[i+1]):

            c = (j+1)*dc

            for k in range(PanelsPerHalfHoop*2): 
                FusePoints.append(Interpolate3D(hooppointlist[i][k],hooppointlist[i+1][k],c))

    dc = 1./PanelsPerSectionlist[-1]
    
    for i in range(PanelsPerSectionlist[-1]-1):

        c = (i+1)*dc
        for j in range(PanelsPerHalfHoop*2): 
            FusePoints.append(Interpolate3D(hooppointlist[-1][j],aftpoint,c))

    FusePoints.append(aftpoint)
    
    FusePoints = np.array(FusePoints)

    panels = []

    for i in range(PanelsPerHalfHoop*2):
        if i < PanelsPerHalfHoop*2-1:
            panels.append([i+2,0,0,i+1])
        else:
            panels.append([1,0,0,i+1])

    totalpanels = sum(PanelsPerSectionlist)*PanelsPerHalfHoop*2

    for i in range(1,sum(PanelsPerSectionlist)-1):
        for j in range(PanelsPerHalfHoop*2):

            if j < PanelsPerHalfHoop*2-1:
                
                firstpoint = (i)*PanelsPerHalfHoop*2 + 2 + j
                secondpoint = (i-1)*PanelsPerHalfHoop*2 +2 + j
                thirdpoint = (i-1)*PanelsPerHalfHoop*2 + 1 + j
                fourthpoint =(i)*PanelsPerHalfHoop*2 + 1 + j
                

            else:
                
                firstpoint = (i)*PanelsPerHalfHoop*2 + 1
                secondpoint = (i-1)*PanelsPerHalfHoop*2 + 1
                thirdpoint = (i-1)*PanelsPerHalfHoop*2 + 1 + j
                fourthpoint = (i)*PanelsPerHalfHoop*2 + j +1
                
            panels.append([firstpoint,secondpoint,thirdpoint,fourthpoint])

    endpoint = len(FusePoints)-1

    for i in range(endpoint-PanelsPerHalfHoop*2,endpoint):
        panels.append([endpoint,i-1,i,endpoint])

    #for k in range(len(panels)):
    #        print panels[k]

    panels = np.array(panels,dtype = int)
              
    return FusePoints,panels

def FuseArea(FusePoints,panels):

    S = 0.0

    for i in range(len(panels)):

        p1x = FusePoints[panels[i][0]][0]
        p1y = FusePoints[panels[i][0]][1]
        p1z = FusePoints[panels[i][0]][2]

        p2x = FusePoints[panels[i][1]][0]
        p2y = FusePoints[panels[i][1]][1]
        p2z = FusePoints[panels[i][1]][2]

        p3x = FusePoints[panels[i][2]][0]
        p3y = FusePoints[panels[i][2]][1]
        p3z = FusePoints[panels[i][2]][2]

        p4x = FusePoints[panels[i][3]][0]
        p4y = FusePoints[panels[i][3]][1]
        p4z = FusePoints[panels[i][3]][2]

        r1 = np.array([p1x-p2x,p1y-p2y,p1z-p2z])
        r2 = np.array([p3x-p2x,p3y-p2y,p3z-p2z])
        r3 = np.array([p3x-p4x,p3y-p4y,p3z-p4z])
        r4 = np.array([p1x-p4x,p1y-p4y,p1z-p4z])

        norm1 = np.cross(r1,r2)
        norm2 = np.cross(r3,r4)

        normfinal = ((norm1+norm2)/2)

        Slocal = np.linalg.norm(normfinal)

        S+= Slocal

    return S

def CreatePlaneGeom(wingarray,Hwingarray,Vwingarray,elipsearray,elipselocarray,PanelsPerHalfHoop,fuselagelength,panelspersection,foilarray,patharray,tfoil):

    coords,uppers,lowers = loadfoils(patharray)

    r2 = wingarray[0]
    r3 = wingarray[1]
    r4 = wingarray[2]
    S = wingarray[3]
    b = wingarray[4]
    xpanels = wingarray[5]
    SpanPanelsPerSection = wingarray[6]
    xwing = wingarray[7]
    fuseelipsewidth = wingarray[8]
    zwing = wingarray[9]
    dihedral = wingarray[10]

    taperh = Hwingarray[0]
    bh = Hwingarray[1]
    Sh = Hwingarray[2]
    xpanelsh = Hwingarray[3]
    SpanPanelsPerSectionh = Hwingarray[4]
    xwingh = Hwingarray[5]
    fuseelipsewidthh = Hwingarray[6]
    zwingh = Hwingarray[7]
    dihedralh = Hwingarray[8]

    taperv = Vwingarray[0]
    bv = Vwingarray[1]
    Sv = Vwingarray[2]
    xpanelsv = Vwingarray[3]
    SpanPanelsPerSectionv = Vwingarray[4]
    xwingv = Vwingarray[5]
    fuseelipsewidthv = Vwingarray[6]
    zwingv = Vwingarray[7]

    
    wingcoords,panels = createwingpoints(r2,r3,r4,b,S,xpanels,SpanPanelsPerSection,xwing,fuseelipsewidth,zwing,dihedral,foilarray,coords,uppers,lowers)
    Htailcoords,Htailpanels = createHtailpoints(taperh,bh,Sh,xpanelsh,SpanPanelsPerSectionh,xwingh,fuseelipsewidthh,zwingh,dihedralh,tfoil,coords,uppers,lowers)
    Vtailcoords,Vtailpanels = createVtailpoints(taperv,bv,Sv,xpanelsv,SpanPanelsPerSectionv,xwingv,fuseelipsewidthv,zwingv,tfoil,coords,uppers,lowers)

    hoopcoords,fusepanels = CreateFuseMesh(elipsearray,elipselocarray,PanelsPerHalfHoop,panelspersection,fuselagelength)

    Sfuse = FuseArea(hoopcoords,fusepanels)

    Vtailcoords[:,1] = 0.0

    planepanels = np.copy(panels)
    planecoords = np.copy(wingcoords)
    
    planepanels = np.append(planepanels,(Htailpanels+len(planecoords)),axis =0)
    planecoords = np.append(planecoords,Htailcoords,axis =0)
    
    planepanels = np.append(planepanels,(Vtailpanels+len(planecoords)),axis =0)
    planecoords = np.append(planecoords,Vtailcoords,axis =0)
    
    #planepanels = np.append(planepanels,(fusepanels+len(planecoords)),axis =0)
    #planecoords = np.append(planecoords,hoopcoords,axis =0)

    return planecoords,planepanels,Sfuse

