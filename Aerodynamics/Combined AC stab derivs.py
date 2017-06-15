###############################################################################
# Aerodynamic center
###############################################################################
# Inputs needed: refx,refy,refz
# need to be defined: wingspan, coord, coords, coordpanels, spanpanels, Force, rho, Vinfx, S 
def calcAC(refx,refy,refz):

    #raw inputs: refx,refy,refz
    def reflength(refx,refy,refz):
        cx = (vortexpoints[:,0,0]-refx)
        cy = (vortexpoints[:,0,1]-refy)+((wingspan/spanpanels)/2)
        cz = (vortexpoints[:,0,2]-refz)
        
        return [cx,cy,cz]
        
    def forcecoeff(Force,rho,Vinfx,S):
        # Cx, X-force coefficient ~ Drag
        Cx = Force[:,0]/(rho*Vinfx*Vinfx*S)
        # Cy, Y-force coefficient ~ Side force
        Cy = Force[:,1]/(rho*Vinfx*Vinfx*S)
        # Cz, Z-force coefficient ~ Lift
        Cz = Force[:,2]/(rho*Vinfx*Vinfx*S)
    
        return [Cx,Cy,Cz]
        
    def momentcoeff(Force,cx,cy,cz,rho,Vinfx,S,wingspan,coord):
        # Cl, rolling moment coefficient
        Mxy = -sum(Force[:,1]*cz)
        Mxz = sum(Force[:,2]*cy)
        l = Mxy+Mxz
        Cl = l/(0.5*rho*Vinfx*Vinfx*S*wingspan)
        # Cm, pitching moment coefficient
        Myx = sum(Force[:,0]*cz)
        Myz = -sum(Force[:,2]*cx)
        m = (Myx+Myz)
        Cm = m/(0.5*rho*Vinfx*Vinfx*S*coord)    
        # Cn, Yawing moment coefficient
        Mzx = -sum(Force[:,0]*cy)
        Mzy = sum(Force[:,1]*cx)
        n = Mzx+Mzy
        Cn = n/(0.5*rho*Vinfx*Vinfx*S*wingspan)
        
        return [Cl,Cm,Cn,l,m,n]

    def momentderiv(Force,cx,cy,cz,rho,Vinfx,S,coord,Cl,Cm,Cn):
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
        Myzmz = -sum(((1.01*Force[:,2])*cx))
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
    
    def calcpanelloc(coords,coordpanels,spanpanels):
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
        
    def actualAC(refx,refy,refz,cx,cy,cz,dclcy,dclcz,dcmcx,dcmcz,dcncx,dcncy,coordlocx,spanlocy,panellocz,coordpanels,spanpanels):
        # Xac, X-location Aerodynamic center
        Xacp = refx+cx*dcmcz+cx*dcncy
        Xac = (sum(Xacp)-(sum(coordlocx)))/spanpanels
        Xacc = Xac/coord
        # Yac, Y-location Aerodynamic center
        Yacp = refy+cy*dclcz+cy*dcncx
        Yac = (sum(Yacp)-(sum(spanlocy)))/coordpanels
        # Zac, Z-location Aerodyanmic center
        Zacp = refz+cz*dclcy+cz*dcmcx
        Zac = (sum(Zacp)-(sum(panellocz)))/spanpanels
            
        return [Xac,Yac,Zac,Xacc]
    
    #General inputs:
        
    reflength = reflength(refx,refy,refz)
    
    cx = reflength[0]
    cy = reflength[1]
    cz = reflength[2]
    
    forcecoeff = forcecoeff(Force,rho,Vinfx,S)
    momentcoeff = momentcoeff(Force,cx,cy,cz,rho,Vinfx,S,wingspan,coord)
    
    Cl = momentcoeff[0]
    Cm = momentcoeff[1]
    Cn = momentcoeff[2]
    l = momentcoeff[3]
    m = momentcoeff[4]
    n = momentcoeff[5]
    
    momentderiv = momentderiv(Force,cx,cy,cz,rho,Vinfx,S,coord,Cl,Cm,Cn)
    
    dclcy = momentderiv[0]
    dclcz = momentderiv[1]
    dcmcx = momentderiv[2]
    dcmcz = momentderiv[3]
    dcncx = momentderiv[4]
    dcncy = momentderiv[5]
    
    panellocs = calcpanelloc(coords,coordpanels,spanpanels)
    
    coordlocx = panellocs[0]
    spanlocy = panellocs[1]
    panellocz = panellocs[2]
    
    AC = actualAC(refx,refy,refz,cx,cy,cz,dclcy,dclcz,dcmcx,dcmcz,dcncx,dcncy,coordlocx,spanlocy,panellocz,coordpanels,spanpanels)
    
    return [AC,l,m,n]
    
def calcstabderivs():#does not use specific inputs
    
    
    def calcstabderv1():
        
        # Run VLM for V+1 here, outputs should be called Forcev1, Vinfxv1    
        
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
        
        return [Cxu,Czu,Cmu]
    
    def calcstabdera1():
        
        # Run VLM for a+1 here, outputs should be called Forcea1, Vinfxa1, FLa1, FDa1
        
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
        
        return [CLa,CDa,Cxa,Cza,Cma]
        
    def calcstabderb1():
        
        # Run VLM for b+1 here, outputs should be called Forceb1, Vinfxb1
        
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
        
        return [Cyb,Clb,Cnb]
    
    def calcstabderad():# need CNha, VhV,deda,ShS and lhc either as inputs or pre-defined
        
        Czad = CNha*VhV*VhV*deda*ShS*lhc
    
        Cmad = -CNha*VhV*VhV*deda*ShS*lhc*lhc
        
        return [Czad,Cmad]
        
    def calcstabderrr():
        dVinfyrr = np.c_[np.zeros((totalpanels, 1)), (-spanlocy+(spanlocy*math.cos(rr))),np.zeros((totalpanels, 1))]
        dVinfzrr =np.c_[np.zeros((totalpanels, 2)), (-spanlocy*math.sin(rr))]
    
        # Run VLM for rr+1 here, outputs should be called Forcerr, Vinfxrr
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
        
        return [Cyp,Clp,Cnp]
    
    def calcstabderpr():
        dVinfxpr = np.c_[(coordlocx-(coordlocx*math.cos(pr))), np.zeros((totalpanels, 2))]        
        dVinfzpr = np.c_[np.zeros((totalpanels, 2)), (coordlocx*math.sin(pr))]
        
    
        # Run VLM for pr+1 here, outputs should be called Forcepr, Vinfxpr
    
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
        
        return [Czq,Cmq]
        
    def calcstabderyr():
        dVinfxyr = np.c_[(spanlocy*math.sin(yr)), np.zeros((totalpanels, 2))]        
        dVinfyyr = np.c_[np.zeros((totalpanels, 1)), (spanlocy-(spanlocy*math.cos(yr))), np.zeros((totalpanels, 1))]
                
        # Run VLM for yr+1 here, outputs should be called Forceyr, Vinfxyr
        
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
        
        return [Clr,Cnr]
    
    AC = calcAC(refx,refy,refz)
    l = AC[1]
    m = AC[2]
    n = AC[3]
    
    Cxu = calcstabderv1()[0]
    Czu = calcstabderv1()[1]
    Cmu = calcstabderv1()[2]
    CLa = calcstabdera1()[0]
    CDa = calcstabdera1()[1]
    Cxa = calcstabdera1()[2]
    Cza = calcstabdera1()[3]
    Cma = calcstabdera1()[4]
    Cyb = calcstabderb1()[0]
    Clb = calcstabderb1()[1]
    Cnb = calcstabderb1()[2]
    Czad = calcstabderad()[0]
    Cmad = calcstabderad()[1]
    Cyp = calcstabderrr()[0]
    Clp = calcstabderrr()[1]
    Cnp = calcstabderrr()[2]
    Czq = calcstabderpr()[0]
    Cmq = calcstabderpr()[1]
    Clr = calcstabderyr()[0]
    Cnr = calcstabderyr()[1]
    
    stabderivs = [Cxu,Czu,Cmu,CLa,CDa,Cxa,Cza,Cma,Cyb,Clb,Cnb,Czad,Cmad,Cyp,Clp,Cnp,Czq,Cmq,Clr,Cnr]    
    
    return [stabderivs]

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

