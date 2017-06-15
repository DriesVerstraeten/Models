## Composite Analysis:
import numpy as np

def moi_ellipse(a,b):
    I_xx = np.pi / 4 * b**2 *(3*a+b)
    I_yy = np.pi / 4 * a**2 *(3*b+a)
    return I_xx,I_yy

def f_x(x,a,b):
    f = (a**2*(np.sin(x)**2)+b**2*(np.cos(x)**2))**0.5
    return f

def arclength(angleboom_top,a,b):
    nodes = np.array([angleboom_top[0],angleboom_top[0]+ \
                 (angleboom_top[1]-angleboom_top[0])/3,\
                  angleboom_top[0]+2./3*(angleboom_top[1]-angleboom_top[0]),\
                  angleboom_top[1]])
    boom_spacing = (angleboom_top[1]-angleboom_top[0])/8* \
                         (f_x(nodes[0],a,b)+ \
                        3*f_x(nodes[1],a,b)+ \
                        3*f_x(nodes[2],a,b)+ \
                        f_x(nodes[3],a,b))
    return boom_spacing

def moi_cockpit(a_ellipse,b_ellipse_top,b_ellipse_bot,d):
    
    angleboom_bot = np.arange(180,360,1)*np.pi/180
    dL_bot = arclength(angleboom_bot,a_ellipse,b_ellipse_bot)    
    
    y_bot = b_ellipse_bot*np.sin(angleboom_bot)
    y_c_bot = b_ellipse_bot+np.sum(y_bot*dL_bot)/(np.shape(y_bot)[0]*dL_bot)
    A_bot = ( a_ellipse + b_ellipse_bot )/2
    
    x_panel = a_ellipse
    y_panel = b_ellipse_bot+d/2
    A_panel = d
    
    angleboom_top = np.arange(0,180,1)*np.pi/180
    dL_top = arclength(angleboom_top,a_ellipse,b_ellipse_top)

    y_top = b_ellipse_top*np.sin(angleboom_top)+d+b_ellipse_bot
    y_c_top = np.sum(y_top*dL_top)/(np.shape(y_top)[0]*dL_top)
    A_top = (a_ellipse+b_ellipse_top)/2

    x_c = 0 
    y_c = (A_bot*y_c_bot+2*A_panel*y_panel+A_top*y_c_top)/(A_bot+A_top+ \
                                                           2*A_panel)
    
    I_xx_top = moi_ellipse(a_ellipse,b_ellipse_top)[0]/2
    I_yy_top = moi_ellipse(a_ellipse,b_ellipse_top)[1]/2
    I_xx_bot = moi_ellipse(a_ellipse,b_ellipse_bot)[0]/2
    I_yy_bot = moi_ellipse(a_ellipse,b_ellipse_bot)[1]/2
    I_xx_panel      = d^3/12
    I_yy_panel      = 0

    I_xx = I_xx_top+A_top*(y_c-y_c_bot)**2+I_xx_bot+A_bot*(y_c-y_c_bot)**2 + \
            2*I_xx_panel+2*A_panel*(y_c-y_panel)**2
    I_yy = I_yy_bot+I_yy_top+I_yy_panel+2*A_panel*(x_c-x_panel)**2
    
    return y_c,I_xx,I_yy

def section_i(a_I,b_I):
    angle = np.arange(361)*np.pi/180
    x = a_I*np.cos(angle)
    y = b_I*np.sin(angle)
    return x,y

def section_ii(a,b_top,b_bot,d,y_c):
    angle_top = np.arange(181)*np.pi/180
    angle_bot = np.arange(180,361)*np.pi/180
    dl  = np.round(arclength(angle_top,a,b_top),1)
    
    x_t = a*np.cos(angle_top)
    y_t = b_top*np.sin(angle_top)+b_bot+d-y_c 
    
    y_d_l   = np.arange(d,0,dl)+b_bot-y_c
    x_d_l   = -a*np.ones(np.shape(y_d)[0])
    
    x = np.concatenate(x_t,x_d)
    y = np.concatenate(y_t,y_d)
    
    x_b = a*np.cos(angle_bot)
    y_b = b_bot*np.sin(angle_bot)
    
    x = np.append(x,x_b)
    y = np.append(y,y_b)
    
    y_d_r = np.arange(0,d,dl)+b_bot-y_c
    x_d_r = -x_d_l
    
    x = np.append(x,x_d_r)
    y = np.append(y,y_d_r)
    return x,y

def bending_i(M_x,M_y,a_i,b_i,material):
    for i in range(np.shape(a_i)[0])
        x = section_i(a_i[i],b_i[i])[0]
        y = section_i(a_i[i],b_i[i])[1]
        I_xx = moi_ellipse(a_i[i],b_i[i])[0]
        I_yy = moi_ellipse(a_i[i],b_i[i])[0]
        sigma_t = M_x/I_xx*y+M_y/I_yy*x 
        
    
    return 