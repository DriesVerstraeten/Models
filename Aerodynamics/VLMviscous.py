import numpy as np
import math
from os import listdir

def CalcWingDrag(planecoords,planepanels,xpanels,data,forces,spanpanels,
               WingStartPanel,V,Rho,foilarray,DynViscosity,fusewidth,b,dihedral):

    ylist = [-b/2,-b/3,-b/6,-fusewidth/2]

    for i in range(len(ylist)):
        ylist[i] = ylist[i]*math.cos(math.radians(dihedral))

    ylist2 = []

    for i in range(len(ylist)):
        ylist2.append(-ylist[-i-1])

    ylist = ylist + ylist2

    totaldrag = 0.0

    for i in range(spanpanels):

        frontpanelindex = i*xpanels+WingStartPanel
        firstpanel = planepanels[frontpanelindex]
        lastpanel = planepanels[frontpanelindex+xpanels-1]

        front1 = planecoords[firstpanel[1]][0]
        front2 = planecoords[firstpanel[2]][0]

        aft1 = planecoords[lastpanel[0]][0]
        aft2 = planecoords[lastpanel[3]][0]

        c1 = aft1-front1
        c2 = aft2-front2

        Clocal = (c1+c2)/2.0

        Reynolds = Rho*V*Clocal/DynViscosity

        panelspan1 = planecoords[firstpanel[1]][1]
        panelspan2 = planecoords[firstpanel[2]][1]

        yloc = (panelspan1+panelspan2)/2

        panelspan = panelspan2-panelspan1

        Totalforce = np.zeros(3)
        
        for j in range(xpanels):
            Totalforce += forces[frontpanelindex+j]

        lift = Totalforce[2]

        localS = panelspan*Clocal

        localCL = lift/(0.5*Rho*V*V*localS)

        k = 0
        while not(not(yloc<ylist[k]) and not(yloc > ylist[k+1])) and (k < len(ylist)+1):
            k +=1

        c = (ylist[k]-yloc)/(ylist[k]-ylist[k+1])

        Cdlocalleft = GetFoilCdfromCL(foilarray[k],data,localCL,Reynolds)
        Cdlocalright = GetFoilCdfromCL(foilarray[k+1],data,localCL,Reynolds)

        cdlocal = Cdlocalleft + c * (Cdlocalright-Cdlocalleft)

        Drag = cdlocal*localS*V*V*Rho*0.5

        totaldrag += Drag

    return totaldrag

def CalcHTailDrag(planecoords,planepanels,xpanels,data,forces,spanpanels,
               WingStartPanel,V,Rho,foil,DynViscosity,fusewidth,bh,dihedralh):

    totaldrag = 0.0

    ylist = [-bh/2,-fusewidth/2]

    for i in range(len(ylist)):
        ylist[i] = ylist[i]*math.cos(math.radians(dihedralh))

    ylist2 = []

    for i in range(len(ylist)):
        ylist2.append(-ylist[-i-1])

    ylist = ylist + ylist2

    for i in range(spanpanels):

        frontpanelindex = i*xpanels+WingStartPanel
        firstpanel = planepanels[frontpanelindex]
        lastpanel = planepanels[frontpanelindex+xpanels-1]

        front1 = planecoords[firstpanel[1]][0]
        front2 = planecoords[firstpanel[2]][0]

        aft1 = planecoords[lastpanel[0]][0]
        aft2 = planecoords[lastpanel[3]][0]

        c1 = aft1-front1
        c2 = aft2-front2

        Clocal = (c1+c2)/2.0

        Reynolds = Rho*V*Clocal/DynViscosity

        panelspan1 = planecoords[firstpanel[1]][1]
        panelspan2 = planecoords[firstpanel[2]][1]

        yloc = (panelspan1+panelspan2)/2

        panelspan = panelspan2-panelspan1

        Totalforce = np.zeros(3)
        
        for j in range(xpanels):
            Totalforce += forces[frontpanelindex+j]

        lift = Totalforce[2]

        localS = panelspan*Clocal

        localCL = lift/(0.5*Rho*V*V*localS)

        k = 0
        while not(not(yloc<ylist[k]) and not(yloc > ylist[k+1])) and (k < len(ylist)+1):
            k +=1

        c = (ylist[k]-yloc)/(ylist[k]-ylist[k+1])

        Cdlocalleft = GetFoilCdfromCL(foil,data,localCL,Reynolds)
        Cdlocalright = GetFoilCdfromCL(foil,data,localCL,Reynolds)

        cdlocal = Cdlocalleft + c * (Cdlocalright-Cdlocalleft)

        Drag = cdlocal*localS*V*V*Rho*0.5

        totaldrag += Drag

    return totaldrag

def CalcVTailDrag(planecoords,planepanels,xpanels,data,forces,spanpanels,
               WingStartPanel,V,Rho,foil,DynViscosity,fusewidth,bv,zwing):

    totaldrag = 0.0

    zlist = [fusewidth/2+zwing,bv+zwing]

    for i in range(spanpanels):

        frontpanelindex = i*xpanels+WingStartPanel
        firstpanel = planepanels[frontpanelindex]
        lastpanel = planepanels[frontpanelindex+xpanels-1]

        front1 = planecoords[firstpanel[1]][0]
        front2 = planecoords[firstpanel[2]][0]

        aft1 = planecoords[lastpanel[0]][0]
        aft2 = planecoords[lastpanel[3]][0]

        c1 = aft1-front1
        c2 = aft2-front2

        Clocal = (c1+c2)/2.0

        Reynolds = Rho*V*Clocal/DynViscosity

        panelspan1 = planecoords[firstpanel[1]][2]
        panelspan2 = planecoords[firstpanel[2]][2]

        zloc = (panelspan1+panelspan2)/2

        panelspan = abs(panelspan2-panelspan1)

        Totalforce = np.zeros(3)
        
        for j in range(xpanels):
            Totalforce += forces[frontpanelindex+j]

        lift = Totalforce[2]

        localS = panelspan*Clocal

        localCL = lift/(0.5*Rho*V*V*localS)

        k = 0
        while not(not(zloc<zlist[k]) and not(zloc > zlist[k+1])) and (k < len(zlist)+1):
            k +=1

        c = (zlist[k]-zloc)/(zlist[k]-zlist[k+1])

        Cdlocalleft = GetFoilCdfromCL(foil,data,localCL,Reynolds)
        Cdlocalright = GetFoilCdfromCL(foil,data,localCL,Reynolds)

        cdlocal = Cdlocalleft + c * (Cdlocalright-Cdlocalleft)

        Drag = cdlocal*localS*V*V*Rho*0.5

        totaldrag += Drag
        
    return totaldrag

def GetFoilCdfromCL(foil,data,CL,Reynolds):
	
    i = 0

    while not(not(Reynolds < data[foil][1][i][0]) and not(Reynolds > data[foil][1][i+1][0])) and (i < (len(data[foil][1])-1)):

        i+=1

    reylow = i
    reyhigh = i+1
    
    j1 = 0
    
    while not(not(CL < data[foil][0][reylow][j1][1]) and not(CL > data[foil][0][reylow][j1+1][1])) and (j1 < len(data[foil][0][reylow])-1):

        j1+=1

    lowjlow = j1
    lowjhigh = j1+1

    j2 = 0
    while not(not(CL < data[foil][0][reyhigh][j2][1]) and not(CL > data[foil][0][reyhigh][j2+1][1])) and (j2 < len(data[foil][0][reylow])-1):

        j2+=1

    highjlow = j2
    highjhigh = j2+1

    lowreyCLlow = data[foil][0][reylow][lowjlow][1]
    lowreyCLhigh = data[foil][0][reylow][lowjhigh][1]

    lowreyCdlow = data[foil][0][reylow][lowjlow][2]
    lowreyCdhigh = data[foil][0][reylow][lowjhigh][2]

    highreyCLlow = data[foil][0][reyhigh][highjlow][1]
    highreyCLhigh = data[foil][0][reyhigh][highjhigh][1]

    highreyCdlow = data[foil][0][reyhigh][highjlow][2]
    highreyCdhigh = data[foil][0][reyhigh][highjhigh][2]

    Cdlow = lowreyCdlow + (CL-lowreyCLlow)*(lowreyCdhigh-lowreyCdlow)/(lowreyCLhigh-lowreyCLlow)
    Cdhigh = highreyCdlow + (CL-highreyCLlow)*(highreyCdhigh-highreyCdlow)/(highreyCLhigh-highreyCLlow)

    lowreyval = data[foil][1][reylow][0]
    highreyval = data[foil][1][reyhigh][0]

    Cd = Cdlow + (Reynolds-lowreyval)*(Cdhigh-Cdlow)/(highreyval-lowreyval)

    return Cd

def GetFuselageDrag(fuselagelength,DynViscosity,S,rho,V):

    Reynolds = fuselagelength*V*rho/DynViscosity

    Cd = 0.075/((math.log(Reynolds,10)-2)**2)

    Drag = Cd*S*rho*V*V*0.5

    return Drag

def SortReyList(reylist):

    for i in range(len(reylist)):
        for j in range(len(reylist)-1):
            
            if(reylist[j][0] > reylist[j+1][0]):
                temp = reylist[j]
                reylist[j] = reylist[j+1]
                reylist[j+1] =temp

    return reylist
    
    
def importpolars(pathlist):

    foils = []

    for i in range(len(pathlist)):

        files = listdir(pathlist[i])

        reynoldslist = []

        for k in range(len(files)):
            f = open(pathlist[i]+ '/' + files[k])
            lines = f.readlines()

            reyline = lines[7]

            reyline = reyline.split("    ")

            a = reyline[2]

            a = a.split(' ')

            if (a[0] == ''):
                del a[0]

            reynolds = float(a[0])*10**float(a[2])

            reynoldslist.append([reynolds,k])

        reylist = SortReyList(reynoldslist)

        foildata = []

        for j in range(len(files)):

            reylistuse = reylist[j][1]
            
            a = np.genfromtxt(pathlist[i]+'/'+files[reylistuse], skip_header = 11)
            foildata.append(a)

        foils.append([foildata,reylist])

    foils = np.array(foils)

    return foils

def ViscousDrag(plane,forces,wingarray,htailarray,vtailarray,data,V,Rho,foilarray,tfoil,DynViscosity,fuseS,alpha,fuselength):

    totalpanels = len(plane[5])
    Force2 = np.empty((totalpanels,3))
    alpharad = math.pi/180.*alpha

    rotmaty = np.matrix([[math.cos(alpharad),0,math.sin(alpharad)],
                         [0,1,0],
                         [-math.sin(alpharad),0,math.cos(alpharad)]])

    for i in range(totalpanels):
        Force2[i] = np.dot(rotmaty,forces[i])

    planepanels = plane[8]
    planecoords = plane[7]

    xpanelswing = wingarray[5]
    xpanelshtail = htailarray[3]
    xpanelsvtail = vtailarray[3]

    spanpanelswing = wingarray[6]*6
    spanpanelhtail = htailarray[4]*2
    spanpanelvtail = vtailarray[4]

    wingstartpanel = plane[6][0]
    htailstartpanel = plane[6][1]
    vtailstartpanel = plane[6][2]

    fusewidth = wingarray[8]
    fusewidthhtail = htailarray[6]
    fusewidthvtail = vtailarray[6]

    b = wingarray[4]
    bh = htailarray[1]
    bv = vtailarray[1]

    dihedral = wingarray[10]
    dihedralh = htailarray[8]

    zVtail = htailarray[7]

    Dragwing = CalcWingDrag(planecoords,planepanels,xpanelswing,data,Force2,spanpanelswing,
               wingstartpanel,V,Rho,foilarray,DynViscosity,fusewidth,b,dihedral)
    DragHtail = CalcHTailDrag(planecoords,planepanels,xpanelshtail,data,Force2,spanpanelhtail,
               htailstartpanel,V,Rho,tfoil,DynViscosity,fusewidthhtail,bh,dihedralh)
    DragVtail = CalcVTailDrag(planecoords,planepanels,xpanelsvtail,data,Force2,spanpanelvtail,
               vtailstartpanel,V,Rho,tfoil,DynViscosity,fusewidthvtail,bv,zVtail)

    Fusedrag = GetFuselageDrag(fuselength,DynViscosity,fuseS,Rho,V)

    totaldrag = Dragwing + DragHtail + DragVtail + Fusedrag

    return totaldrag

