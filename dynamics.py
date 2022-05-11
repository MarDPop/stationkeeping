# -*- coding: utf-8 -*-
import numpy as np
import solarsystem

def dxdt(x,t,mu,mu1):
    dx1 = x[0] + mu
    dx2 = x[0] - mu1
    d = x[1]*x[1] + x[2]*x[2]
    r1_sq = dx1*dx1 + d
    r2_sq = dx2*dx2 + d
    
    r1 = np.sqrt(r1_sq)
    r2 = np.sqrt(r2_sq)
    
    g1 = mu1/(r1*r1_sq)
    g2 = mu/(r2*r2_sq)
    g = g1 + g2
    
    x_dd = x[0] + 2*x[4] - dx1*g1 - dx2*g2
    y_dd = x[1]*(1 - g) - 2*x[3]
    z_dd = -x[2]*g
    
    return np.array([x[3], x[4], x[5], x_dd, y_dd,z_dd])

def dxdt_earth_moon_sun(x,t,jd0): 
    jd = jd0 + t/86400.0
    #EARTH
    r = np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    r3 = r*r*r
    x_dd = solarsystem.EARTH_MU*x[0]/r3
    y_dd = solarsystem.EARTH_MU*x[1]/r3
    z_dd = solarsystem.EARTH_MU*x[2]/r3
    #MOON
    eci = solarsystem.getMoonECI(jd)
    dx = [eci[0] - x[0],eci[1] - x[1],eci[2] - x[2]]
    r = np.sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])
    r3 = r*r*r
    x_dd += solarsystem.MOON_MU*dx[0]/r3
    y_dd += solarsystem.MOON_MU*dx[1]/r3
    z_dd += solarsystem.MOON_MU*dx[2]/r3
    #SUN
    eci = solarsystem.getSunECI(jd)
    dx = [eci[0] - x[0],eci[1] - x[1],eci[2] - x[2]]
    r = np.sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])
    r3 = r*r*r
    x_dd += solarsystem.SUN_MU*dx[0]/r3
    y_dd += solarsystem.SUN_MU*dx[1]/r3
    z_dd += solarsystem.SUN_MU*dx[2]/r3
    
    return np.array([x[3], x[4], x[5], x_dd, y_dd,z_dd])

def dxdt_earth_moon_sun_solarpressure(x,t,jd0,bc): 
    jd = jd0 + t/86400.0
    #EARTH
    r = np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    r3 = r*r*r
    x_dd = -solarsystem.EARTH_MU*x[0]/r3
    y_dd = -solarsystem.EARTH_MU*x[1]/r3
    z_dd = -solarsystem.EARTH_MU*x[2]/r3
    #MOON
    eci = solarsystem.getMoonECI(jd)
    dx = [x[0] - eci[0],x[1] - eci[1], x[2] - eci[2]]
    r = np.sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])
    r3 = r*r*r
    x_dd -= solarsystem.MOON_MU*dx[0]/r3
    y_dd -= solarsystem.MOON_MU*dx[1]/r3
    z_dd -= solarsystem.MOON_MU*dx[2]/r3
    #SUN
    eci = solarsystem.getSunECI(jd)
    dx = [x[0] - eci[0],x[1] - eci[1], x[2] - eci[2]]
    r = np.sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])
    r3 = r*r*r
    x_dd -= solarsystem.SUN_MU*dx[0]/r3
    y_dd -= solarsystem.SUN_MU*dx[1]/r3
    z_dd -= solarsystem.SUN_MU*dx[2]/r3
    
    a = solarsystem.SOLAR_PRESSURE_AU*(solarsystem.AU/r)**2/bc
    x_dd += a*dx[0]/r
    y_dd += a*dx[1]/r
    z_dd += a*dx[2]/r
    
    return np.array([x[3], x[4], x[5], x_dd, y_dd,z_dd])

def getA(x,mu,mu1):
    A = np.zeros((6,6))
    A[0][3] = 1
    A[1][4] = 1
    A[2][5] = 1
    A[3][4] = 2
    A[4][3] = -2
    
    dx1 = x[0] + mu
    dx2 = x[0] - mu1
    d = x[1]**2 + x[2]**2
    r1_sq = dx1*dx1 + d
    r2_sq = dx2*dx2 + d
    r1 = np.sqrt(r1_sq)
    r2 = np.sqrt(r2_sq)
    
    c1 = mu1/(r1*r1_sq)
    c2 = mu/(r2*r2_sq)
    
    a = c1 + c2
    
    c1 *= 3/r1_sq
    c2 *= 3/r2_sq
    
    b = c1 + c2
    
    c1 *= dx1
    c2 *= dx2

    d = c1 + c2
    
    yz = x[1]*x[2]
    
    A[3][0] = c1*dx1 + c2*dx2 - a + 1
    A[3][1] = x[1]*d
    A[3][2] = x[2]*d
    
    A[4][0] = x[1]*d
    A[4][1] = x[1]*x[1]*b - a + 1
    A[4][2] = yz*b
    
    A[5][0] = x[2]*d
    A[5][1] = yz*b
    A[5][2] = x[2]*x[2]*b - a

    return A

def getA_finite_diff(f,x0,t,dx):
    A = np.zeros((6,6))
    
    h0 = f(x0,t)
    
    for i in range(6):
        x = x0.copy()
        x[i] += dx[i]
        h = f(x,t)
        dh = (h - h0)/dx[i]
        for j in range(6):
            A[j,i] = dh[j]
            
    return A
 
def getCL1(yL,mu,n):
    n1 = n+1
    return 1/yL**3* (mu + (-1)**n*(1 - mu)*yL**n1/(1 - yL)**n1)

def getCL2(yL,mu,n):
    n1 = n+1
    return 1/yL**3* ((-1)**n*mu + (-1)**n*(1 - mu)*yL**n1/(1 + yL)**n1)

def getL1(mu):
    mu1 = 1 - mu
    L1 = 1 - np.cbrt(mu/3)
    A = -1
    B = 1
    
    for i in range(10):
        f = L1 + A*mu1/(L1 + mu)**2 + B*mu/(L1 - mu1)**2
        df = 1 - 2*A*mu1/(L1 + mu)**3 - 2*B*mu/(L1 - mu1)**3
        L1 = L1 - 0.5*f/df
        
    return L1

def getInitialStateFirstOrder(xL1,mu):

    return 0
    

def getInitialState(xL1,mu,mean_motion,sma,Az0,time = 0):
    # https://commons.erau.edu/cgi/viewcontent.cgi?article=1565&context=edt
    # https://hal.archives-ouvertes.fr/hal-00312910v2/document
    # https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1980CeMec..22..241R&defaultprint=YES&filetype=.pdf
    mu1 = 1 - mu
    yL = xL1 - mu1 # should be negative for L1
    yLA = abs(yL)
    c2 = mu/yLA**3 + mu1/(1 - yLA)**3
    c3 = getCL1(yLA,mu,3)
    c4 = getCL1(yLA,mu,4)
    tmp = np.sqrt(c2*(9*c2 - 8))
    gamma2 = (2 - c2 + tmp)/2
    gamma = np.sqrt(gamma2)
    
    k = 2*gamma/(1 + gamma2 - c2)
    k2 = k*k
    
    d1 = 3*gamma2/k*(k*(6*gamma2-1) - 2*gamma)
    d2 = 8*gamma2/k*(k*(11*gamma2-1) - 2*gamma)
    
    a21 = 0.75*c3*(k2 - 2)/(1 + 2*c2)
    a22 = 0.75*c3/(1 + 2*c2)
    a23 = -0.75*c3*gamma/(k*d1)*(3*k**3*gamma - 6*k*(k-gamma) + 4)
    a24 = -0.75*c3*gamma/(k*d1)*(2 + 3*k*gamma)
    b21 = -1.5*c3*gamma/d1*(3*k*gamma-4)
    b22 = 3*c3*gamma/d1
    d21 = -c3/(2*gamma2)
    
    tmp = 3/(64*gamma2)
    d31 = tmp*(4*c3*a24 + c4)
    d32 = tmp*(4*c3*(a23 - d21) + c4*(4 + k2))
    
    tmp = (9*gamma2 + 1 - c2)/2
    a31 = -2.25*gamma/d2*(4*c3*(k*a23 - b21) + k*c4*(4 + k2)) + tmp/d2*(3*c3*(2*a23 - k*b21) + c4*(2 + 3*k2))
    a32 = -1/d2*(2.25*gamma*(4*c3*(k*a24 - b22) + k*c3) + 3*tmp*(c3*(k*b22 + d21 - 2*a24) - c4))
    
    tmp = 0.375*(9*gamma2 + 1 + c2)
    b31 = (3*gamma*(3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k2)) + tmp*(4*c3*(k*a23 - b21) + k*c4*(4 + k2)))/d2
    b32 = (9*gamma*(c3*(k*b22 + d21 - 2*a24) - c4) + tmp*(4*c3*(k*a24 - b22) + k*c4))/d2
    
    a1 = -1.5*c3*(2*a21 + a23 + 5*d21) - 0.375*c4*(12 - k2)
    a2 = 1.5*c3*(a24 - 2*a22) + 1.125*c4
    
    tmp = 1/(2*gamma*(gamma*(1 + k2) - 2*k))
    s1 = tmp*(1.5*c3*(2*a21*(k2 - 2) - a23*(k2 + 2) - 2*k*b21) - 0.375*c4*(k2*(3*k2 - 8) + 8) )
    s2 = tmp*(1.5*c3*(2*a22*(k2 - 2) + a24*(k2 + 2) + 2*k*b22 + 5*d21) + 0.375*c4*(12 - k2))
    
    l1 = a1 + 2*gamma2*s1
    l2 = a2 + 2*gamma2*s2
    
    delta = gamma2 - c2
    
    a1 = 1/yLA
    r1 = ((1-mu) - xL1)*sma
    scale = yLA
    Az = Az0/r1
    m = 1 # or 3
    dn = 2 - m
    
    n1 = np.sqrt(yLA**3) # normalized
    s = n1*time
    
    Az2 = Az*Az
    Ax2 = -(Az2*l2 + delta)/l1
    Ax = np.sqrt(Ax2)
    
    omega2 = s1*Ax2 + s2*Az2
    omega = 1 + omega2
    tau = omega*s
    
    tau1 = gamma*tau
    
    x = a21*Ax2 + a22*Az2 - Ax*np.cos(tau1) + (a23*Ax2 - a24*Az2)*np.cos(2*tau1) + Ax*(a31*Ax2 - a32*Az2)*np.cos(3*tau1)  
    z = dn* ( Az*np.cos(tau1) + d21*Ax*Az*(np.cos(2*tau1) - 3) + Az*(d32*Ax2 - d31*Az2)*np.cos(3*tau1) )
    
    v = gamma*omega*n1*(k*Ax*np.cos(tau) + 2*(b21*Ax2 - b22*Az2)*np.cos(2*tau) + 3*Ax*(b31*Ax2 - b32*Az2)*np.cos(3*tau) )
     
    ydot = xL1 - (xL1*n1 - scale*v)/n1
    x0 = np.array([x*scale + xL1,0,z*scale,0,ydot,0])
    return x0

def getInitialStateWithPhase(xL1,mu,mean_motion,sma,Az0,phi,time = 0):
    # https://commons.erau.edu/cgi/viewcontent.cgi?article=1565&context=edt
    # https://hal.archives-ouvertes.fr/hal-00312910v2/document
    # https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1980CeMec..22..241R&defaultprint=YES&filetype=.pdf
    mu1 = 1 - mu
    yL = xL1 - mu1 # should be negative for L1
    yLA = abs(yL)
    c2 = mu/yLA**3 + mu1/(1 - yLA)**3
    c3 = getCL1(yLA,mu,3)
    c4 = getCL1(yLA,mu,4)
    tmp = np.sqrt(c2*(9*c2 - 8))
    gamma2 = (2 - c2 + tmp)/2
    gamma = np.sqrt(gamma2)
    
    omegaP = np.sqrt((2 - c2 + tmp)/2)
    omegaV = np.sqrt(c2)
    k = 2*gamma/(1 + gamma2 - c2)
    k2 = k*k
    
    d1 = 3*gamma2/k*(k*(6*gamma2-1) - 2*gamma)
    d2 = 8*gamma2/k*(k*(11*gamma2-1) - 2*gamma)
    
    a21 = 0.75*c3*(k2 - 2)/(1 + 2*c2)
    a22 = 0.75*c3/(1 + 2*c2)
    a23 = -0.75*c3*gamma/(k*d1)*(3*k**3*gamma - 6*k*(k-gamma) + 4)
    a24 = -0.75*c3*gamma/(k*d1)*(2 + 3*k*gamma)
    b21 = -1.5*c3*gamma/d1*(3*k*gamma-4)
    b22 = 3*c3*gamma/d1
    d21 = -c3/(2*gamma2)
    
    tmp = 3/(64*gamma2)
    d31 = tmp*(4*c3*a24 + c4)
    d32 = tmp*(4*c3*(a23 - d21) + c4*(4 + k2))
    
    tmp = (9*gamma2 + 1 - c2)/2
    a31 = -2.25*gamma/d2*(4*c3*(k*a23 - b21) + k*c4*(4 + k2)) + tmp/d2*(3*c3*(2*a23 - k*b21) + c4*(2 + 3*k2))
    a32 = -1/d2*(2.25*gamma*(4*c3*(k*a24 - b22) + k*c3) + 3*tmp*(c3*(k*b22 + d21 - 2*a24) - c4))
    
    tmp = 0.375*(9*gamma2 + 1 + c2)
    b31 = (3*gamma*(3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k2)) + tmp*(4*c3*(k*a23 - b21) + k*c4*(4 + k2)))/d2
    b32 = (9*gamma*(c3*(k*b22 + d21 - 2*a24) - c4) + tmp*(4*c3*(k*a24 - b22) + k*c4))/d2
    
    a1 = -1.5*c3*(2*a21 + a23 + 5*d21) - 0.375*c4*(12 - k2)
    a2 = 1.5*c3*(a24 - 2*a22) + 1.125*c4
    
    tmp = 1/(2*gamma*(gamma*(1 + k2) - 2*k))
    s1 = tmp*(1.5*c3*(2*a21*(k2 - 2) - a23*(k2 + 2) - 2*k*b21) - 0.375*c4*(k2*(3*k2 - 8) + 8) )
    s2 = tmp*(1.5*c3*(2*a22*(k2 - 2) + a24*(k2 + 2) + 2*k*b22 + 5*d21) + 0.375*c4*(12 - k2))
    
    l1 = a1 + 2*gamma2*s1
    l2 = a2 + 2*gamma2*s2
    
    delta = gamma2 - c2
    
    a1 = 1/yLA
    r1 = ((1-mu) - xL1)*sma
    scale = yLA
    Az = Az0/r1
    m = 1 # or 3
    zeta = phi + m*np.pi/2
    dn = 2 - m
    
    n1 = np.sqrt(yLA**3) # normalized
    s = n1*time
    
    Az2 = Az*Az
    Ax2 = -(Az2*l2 + delta)/l1
    Ax = np.sqrt(Ax2)
    
    omega2 = s1*Ax2 + s2*Az2
    omega = 1 + omega2
    tau = omega*s
    
    tau1 = gamma*tau + phi
    
    x = a21*Ax2 + a22*Az2 - Ax*np.cos(tau1) + (a23*Ax2 - a24*Az2)*np.cos(2*tau1) + Ax*(a31*Ax2 - a32*Az2)*np.cos(3*tau1)  
    y = k*Ax*np.sin(tau1) + (b21*Ax2 - b22*Az2)*np.sin(2*tau1) + Ax*(b31*Ax2 - b32*Az2)*np.sin(3*tau1)
    z = dn* ( Az*np.cos(tau1) + d21*Ax*Az*(np.cos(2*tau1) - 3) + Az*(d32*Ax2 - d31*Az2)*np.cos(3*tau1) )
    
    v = gamma*omega*n1*(k*Ax*np.cos(tau) + 2*(b21*Ax2 - b22*Az2)*np.cos(2*tau) + 3*Ax*(b31*Ax2 - b32*Az2)*np.cos(3*tau) )
     
    ydot = xL1 - (xL1*n1 - scale*v)/n1
    xdot = 0
    zdot = 0
    x0 = np.array([x*scale + xL1,y*scale,z*scale,xdot,ydot,zdot])
    return x0


def cr3bp2eci(x,jd):
    moon = solarsystem.getMoonECI(jd)
    pass