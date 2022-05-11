import numpy as np
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import dynamics
import solarsystem

def plots(x0,t_final,sma):
    mu = 0.012155650403206974
    mu1 = 0.987844349596793
    mean_motion = 2.6590930417337446e-06
    L1 = dynamics.getL1(mu)
    
    print("running to time (s): " + str(t_final))
    x = x0.copy()
    
    t = 0
    t_record = 0
    record_interval = 0.02
    dt = 0.0001
    t_final = t_final * mean_motion

    # RK Constants
    dt2 = dt/2
    dt3 = dt/3
    dt6 = dt/6
    dt8 = dt/8

    ts = []
    X = []
    Y = []
    Z = []
    while(t < t_final):
        xn = x.copy()
        
        if t > t_record:
            X.append(xn[0]*sma)
            Y.append(xn[1]*sma)
            Z.append(xn[2]*sma)
            ts.append(t)
            t_record += record_interval
            
        k0 = dynamics.dxdt(x,t,mu,mu1)
            
        x = xn + k0*dt3
        
        k1 = dynamics.dxdt(x,t + dt3,mu,mu1)
        
        x = xn + (k0 + k1)*dt6

        k2 = dynamics.dxdt(x,t + dt3,mu,mu1)*3

        x = xn + (k0 + k2)*dt8

        k3 = dynamics.dxdt(x,t + dt2,mu,mu1)*4

        x = xn + (k0 - k2 + k3)*dt2

        k4 = dynamics.dxdt(x,t + dt,mu,mu1)

        x = xn + (k0 + k3 + k4)*dt6
        
        t  += dt

    print('plotting')
       
    RE = 6371.1
    RM = 1731.1
    
    fig, ax = plt.subplots(2)
    
    circle1 = plt.Circle((-mu*sma, 0), RE, color='b')
    circle2 = plt.Circle(((1-mu)*sma, 0), RM, color='r')
    ax[0].add_patch(circle1)
    ax[0].add_patch(circle2)
    ax[0].scatter(L1*sma,0,marker='*',color="y")

    circle1 = plt.Circle((-mu*sma, 0), RE, color='b')
    circle2 = plt.Circle(((1-mu)*sma, 0), RM, color='r')
    ax[1].add_patch(circle1)
    ax[1].add_patch(circle2)
    ax[1].scatter(L1*sma,0,marker='*',color="y")

    ax[0].plot(X,Y)
    ax[0].set_aspect('equal', 'box')

    ax[1].plot(X,Z)
    ax[1].set_aspect('equal', 'box')

    fig = plt.figure()
    ax3d = plt.axes(projection='3d')

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    xS = np.cos(u)*np.sin(v)
    yS = np.sin(u)*np.sin(v)
    zS = np.cos(v)

    xE = RE*xS - mu*sma
    xM = RM*xS + mu1*sma
    ax3d.plot_wireframe(xE, RE*yS, RE*zS, color="b")
    ax3d.plot_wireframe(xM, RM*yS, RM*zS, color="r")
    ax3d.scatter(L1*sma,0,0,marker='*',color="y")

    ax3d.plot(X,Y,Z)
    ax3d.set_xlim([-0.2*sma,1*sma])
    ax3d.set_ylim([-0.6*sma,0.6*sma])
    ax3d.set_zlim([-0.6*sma,0.6*sma])

    plt.show()
    
def plotstationkeeping(x_station,t_station,x_nominal,dVs,dVtime,sma):
    mu = 0.012155650403206974
    mu1 = 0.987844349596793
    # mean_motion = np.sqrt((solarsystem.MOON_MU + solarsystem.EARTH_MU)/(sma**3))
    L1 = dynamics.getL1(mu)

    print('plotting')
       
    RE = 6371.1/sma
    RM = 1731.1/sma

    fig = plt.figure()
    ax = plt.axes(projection='3d')

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    xS = np.cos(u)*np.sin(v)
    yS = np.sin(u)*np.sin(v)
    zS = np.cos(v)

    xE = RE*xS - mu
    xM = RM*xS + mu1
    ax.plot_wireframe(xE, RE*yS, RE*zS, color="b")
    ax.plot_wireframe(xM, RM*yS, RM*zS, color="r")
    ax.scatter(L1,0,0,marker='*',color="y")

    X = []
    Y = []
    Z = []
    for i in range(len(x_nominal)):
        X.append(x_nominal[i][0])
        Y.append(x_nominal[i][1])
        Z.append(x_nominal[i][2])
    ax.plot(X,Y,Z,'b')
    X = []
    Y = []
    Z = []
    for i in range(len(x_station)):
        X.append(x_station[i][0])
        Y.append(x_station[i][1])
        Z.append(x_station[i][2])
    ax.plot(X,Y,Z,'k')
    
    tidx = 0
    for i in range(len(dVtime)):
        t = dVtime[i]
        while t_station[tidx] < t:
            tidx += 1
            
        ax.scatter(x_station[tidx][0],x_station[tidx][1],x_station[tidx][2],marker='o',color="r")

    plt.show()

def patch(ax, x, y, z, color):
    pc = Poly3DCollection([list(zip(x,y,z))])       # Create PolyCollection from coords
    pc.set_facecolor(color)                         # Set facecolor to mapped value
    pc.set_edgecolor('k')                           # Set edgecolor to black
    ax.add_collection3d(pc)                         # Add PolyCollection to axes
    return pc

def dot3(u,v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def heatmap(T, X, Y, Z, sma):
    mu1 = 0.987844349596793
    RM = 1731.1/sma

    nLon = 72
    nLat = 36

    dLon = 2*np.pi/nLon
    dLat = np.pi/nLat

    pointCenter = []
    pointPatch1 = []
    pointPatch2 = []
    for i in range(nLon):
        lonStart = i*dLon
        lonEnd = (i+1)*dLon
        lon = (lonStart + lonEnd)*0.5
        for j in range(nLat):
            latStart = j*dLat
            latEnd = (j+1)*dLat
            lat = (latStart + latEnd)*0.5

            xS = RM*np.cos(lon)*np.sin(lat) + mu1
            yS = RM*np.sin(lon)*np.sin(lat)
            zS = RM*np.cos(lat)
            
            pointCenter.append((xS,yS,zS))

            x1 = np.cos(lonStart)*np.sin(latStart)
            y1 = np.sin(lonStart)*np.sin(latStart)
            z1 = np.cos(latStart)

            x2 = np.cos(lonEnd)*np.sin(latStart)
            y2 = np.sin(lonEnd)*np.sin(latStart)
            z2 = np.cos(latStart)

            x3 = np.cos(lonEnd)*np.sin(latEnd)
            y3 = np.sin(lonEnd)*np.sin(latEnd)
            z3 = np.cos(latEnd)

            x4 = np.cos(lonStart)*np.sin(latEnd)
            y4 = np.sin(lonStart)*np.sin(latEnd)
            z4 = np.cos(latEnd)
            
            pointPatch1.append([[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]])
            pointPatch2.append([[x4,x1,x2],[y4,y1,y2],[z4,z1,z2]])

    tIdx = range(len(T))
    pIdx = range(len(pointCenter))
    pointCount = np.zeros((len(pointCenter),1))
    for i in tIdx:
        x = X[i]
        y = Y[i]
        z = Z[i]
        for j in pIdx:
            p = pointCenter[j]
            dot = p[0]*(p[0] - x) + p[1]*(p[1] - y) + p[2]*(p[2] - z)

            if dot < 0:
                pointCount[j] += 1


    maxCount = 0
    for j in pIdx:
        if pointCount[j] > maxCount:
            maxCount = pointCount[j]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for j in pIdx:
        cVal = pointCount[j]*1.0/maxCount
        c = (cVal,cVal,cVal)
        patch(ax, pointPatch1[j][0], pointPatch1[j][1], pointPatch1[j][2], c)

    plt.show()
    
    
