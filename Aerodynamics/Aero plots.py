import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import time

def calclocs(coordlocx,spanlocy,coordpanels,spanpanels,totalpanels):
    rspanlocy = (np.add.reduceat(spanlocy, np.arange(0, len(spanlocy), (totalpanels/spanpanels))))/coordpanels
    rcoordlocx = coordlocx[:10]

    return [rcoordlocx,rspanocy]
    
def roundup10(x):
    return int(math.ceil(x / 10.0)) * 10

def roundup1(x):
    return int(math.ceil(x / 1.0)) * 1

def plot3dlift(dFz):
    # initiate plot   
    fig = plt.figure(figsize=plt.figaspect(0.5)*1.25)
    ax = Axes3D(fig)

    # Make data.
    X = rcoordlocx
    Y = rspanlocy
    X, Y = np.meshgrid(X, Y)
    Fzplot = dFz
    dFzplot = []
    for i in range(rspanlocy.size):
        dFzplot.append(Fzplot[i*10:i*10+10])
    Z = np.array(dFzplot)
    Zw = np.zeros((rspanlocy.size,rcoordlocx.size))
    color_dimension = Z # change to desired fourth dimension
    minn, maxx = -np.amax(Z), np.amax(Z)
    norm = matplotlib.colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors = fcolors,
            linewidth=0.1, vmin=minn, vmax=maxx, antialiased=True, shade=False)
    surf.set_edgecolor('k')
    surf2 = ax.plot_surface(X, Y, Zw, rstride=1, cstride=1, color='black', 
                          alpha = 0.33, linewidth=0, antialiased=True)
    # Customize the axes.
    ax.set_zlim(-roundup10(np.amax(abs(Z))), roundup10(np.amax(abs(Z))))
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_xlabel('coordloc [m]')
    ax.set_ylabel('spanloc [m]')
    ax.set_zlabel('delta Lift [N]')
        
    # plot figure and set view parameters
    plt.colorbar(m)
    ax.view_init(15,-60)
    
def plot3ddrag(dFz):
    # initiate plot   
    fig = plt.figure(figsize=plt.figaspect(0.5)*1.25)
    ax = Axes3D(fig)

    # Make data.
    X = rcoordlocx
    Y = rspanlocy
    X, Y = np.meshgrid(X, Y)
    Fzplot = dFz
    dFzplot = []
    for i in range(rspanlocy.size):
        dFzplot.append(Fzplot[i*10:i*10+10])
    Z = np.array(dFzplot)
    Zw = np.zeros((rspanlocy.size,rcoordlocx.size))
    color_dimension = Z # change to desired fourth dimension
    minn, maxx = np.amin(Z), -np.amin(Z)
    norm = matplotlib.colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)

    # Plot the surface.
    surf = ax.plot_surface(Z, Y, Zw, rstride=1, cstride=1, facecolors = fcolors,
            linewidth=0.1, vmin=minn, vmax=maxx, antialiased=True, shade=False)
    surf.set_edgecolor('k')
    surf2 = ax.plot_surface(X, Y, Zw, rstride=1, cstride=1, color='black', 
                           alpha = 0.33, linewidth=0, antialiased=True)
    # Customize the axes.
    ax.set_zlim(-2,2)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_xlabel('delta Drag [m]')
    ax.set_ylabel('spanloc [m]')
    ax.set_zlabel('heightloc [N]')
        
    # plot figure and set view parameters
    plt.colorbar(m)
    ax.view_init(25,-60)
    
def plotalpha(alphamin,alphamax,da):
    alphalist = []
    CLlist = []
    CDlist = []
    Cmlist = []
    palpha = alphamin
    
    while palpha < alphamax :
        
        # VLM loop for CL and CD ##################################################    
        pVinfz = V*math.sin(math.pi/180.*palpha)
        pVinfy = 0
        pVinfx = V*math.cos(math.pi/180.*palpha)
    
        pVvec = np.full((totalpanels,3),np.array([pVinfx,pVinfy,pVinfz]))
    
        pb = np.empty(totalpanels)
        for i in range(totalpanels):
            pb[i] = np.dot(pVvec[i],normal[i])
    
        pA = np.empty((totalpanels,totalpanels))
    
        CalcInflMatHorseShoe(pA,vortexpoints,controlpoints)
        calcVind(Vindmat,vortexpoints)
        px = np.linalg.solve(pA,-pb)
    
        pVinducedarr = Vinduced(Vindmat,px)
    
        pVinducedarr2 = pVvec + pVinducedarr
    
        pForce = ((np.cross(pVinducedarr2,(vortexpoints[:,1]-vortexpoints[:,0]))*rho).T*px).T
        pForce2 = np.empty((totalpanels,3))
        palpharad = math.pi/180.*palpha
    
        protmaty = np.matrix([[math.cos(palpharad),0,math.sin(palpharad)],
                         [0,1,0],
                         [-math.sin(palpharad),0,math.cos(palpharad)]])
    
        for i in range(totalpanels):
            pForce2[i] = np.dot(protmaty,pForce[i])
    
        pFx = sum(pForce[:,0])
        pFy = sum(pForce[:,1])
        pFz = sum(pForce[:,2])
    
        pFL = sum(pForce2[:,2])
        pFD = sum(pForce2[:,0])
    
        pS = coord*wingspan
        pA = wingspan*wingspan/S
    
        pCL = pFL/(0.5*rho*V*V*S)
        pCD = pFD/(0.5*rho*V*V*S)
    
        pe = (pCL*pCL)/(math.pi*pA*pCD)
        
        # AC loop for Cm     
        AC = calcAC(refx,refy,refz)
               
        CLlist.append([pCL])
        CDlist.append([pCD])
        Cmlist.append([AC[5]])
        alphalist.append([palpha])
        palpha = palpha + da

    CLarray = np.array(CLlist, dtype=np.float)
    CDarray = np.array(CDlist, dtype=np.float)
    
    LDlist = CLarray/CDarray


start = time.clock()

plotalpha(-5,15,2)

# CL-alpha
plt.figure(1)    
plt.plot(alphalist,CLlist)
plt.axis([alphamin,alphamax,-0.5,2])
plt.xlabel('alpha')
plt.ylabel('CL')
plt.grid(True)
plt.savefig("CL-alpha.png")

# CD-alpha
plt.figure(2)    
plt.plot(alphalist,CDlist)
plt.axis([alphamin,alphamax,0,0.05])
plt.xlabel('alpha')
plt.ylabel('CD')
plt.grid(True)
plt.savefig("CD-alpha.png")

# CL-CD
plt.figure(3)
plt.plot(CDlist,CLlist)
plt.axis([-0.01,0.1,-1,2])
plt.xlabel('CD')
plt.ylabel('CL')
plt.grid(True)
plt.savefig("CL-CD.png")

# Cm-alpha
plt.figure(4)
plt.plot(alphalist,Cmlist)
plt.axis([alphamin,alphamax,-0.5,0.5])
plt.xlabel('alpha')
plt.ylabel('Cm')
plt.grid(True)
plt.savefig("Cm-alpha.png")

# L/D-alpha
plt.figure(5)
plt.plot(alphalist,LDlist)
plt.axis([alphamin,alphamax,-150,150])
plt.xlabel('alpha')
plt.ylabel('L/D')
plt.grid(True)
plt.savefig("LD-alpha.png")

# Lift distribution
plot3dlift(Force[:,2])
plt.savefig("lift_distribution.png")

# Drag distribution
plot3ddrag(Force[:,0])
plt.savefig("drag_distribution.png")

# delta lift distribution v+1
plot3dlift(dFzv1)
plt.savefig("deltalift_distribution_v+1.png")

# delta lift distribution a+1
plot3dlift(dFza1)
plt.savefig("deltalift_distribution_a+1.png")

# delta lift distribution b+1
plot3dlift(dFzb1)
plt.savefig("deltalift_distribution_b+1.png")

# delta lift distribution rr
plot3dlift(dFzrr)
plt.savefig("deltalift_distribution_rr.png")

# delta lift distribution pr
plot3dlift(dFzpr)
plt.savefig("deltalift_distribution_pr.png")

# delta lift distribution yr
plot3dlift(dFzyr)
plt.savefig("deltalift_distribution_yr.png")

# delta drag distribution v+1
plot3ddrag(dFxv1)
plt.savefig("deltadrag_distribution_v+1.png")

# delta drag distribution a+1
plot3ddrag(dFxa1)
plt.savefig("deltadrag_distribution_a+1.png")

# delta drag sitribution b+1
plot3ddrag(dFxb1)
plt.savefig("deltadrag_distribution_b+1.png")

# delta drag distribution rr
plot3ddrag(dFxrr)
plt.savefig("deltadrag_distribution_rr.png")

# delta drag distribution pr
plot3ddrag(dFxpr)
plt.savefig("deltadrag_distribution_pr.png")

# delta drag distribution yr
plot3ddrag(dFxyr)
plt.savefig("deltadrag_distribution_yr.png")

end = time.clock()

print end-start