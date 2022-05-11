import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import dynamics

mu = 0.01215058426954
mu1 = 1-mu
mean_motion = 2.6590930417337446e-06
L1 = dynamics.getL1(mu)

def getX(t,ts,x,y,z):
    lo = 0
    hi = len(ts) - 1;
    mid = int((lo + hi)/2)
    while(lo != mid):
        if t > ts[mid]:
            lo = mid
        else:
            hi = mid
        mid = int((lo + hi)/2)
        
    dt = (t - ts[mid])/(ts[mid+1] - ts[mid])
    
    return np.array([x[mid] + (x[mid+1] - x[mid])*dt,y[mid] + (y[mid+1] - y[mid])*dt,z[mid] + (z[mid+1] - z[mid])*dt])
    

f = open("test_orbit2", "r")
data = f.readlines()

X = []
Y = []
Z = []
X_EML1 = []
Y_EML1 = []
Z_EML1 = []
T = []
for line in data:
    row = line.split()
    T.append(float(row[0]))
    X.append(float(row[1]))
    Y.append(float(row[2]))
    Z.append(float(row[3]))
    X_EML1.append(float(row[4]))
    Y_EML1.append(float(row[5]))
    Z_EML1.append(float(row[6]))
    
f.close()

f = open("test_earth", "r")
data = f.readlines()

X_E = []
Y_E = []
Z_E = []
T_E = []
for line in data:
    row = line.split()
    T_E.append(float(row[0]))
    X_E.append(float(row[1]))
    Y_E.append(float(row[2]))
    Z_E.append(float(row[3]))
f.close()

f = open("test_moon", "r")
data = f.readlines()

X_M = []
Y_M = []
Z_M = []
T_M = []
for line in data:
    row = line.split()
    T_M.append(float(row[0]))
    X_M.append(float(row[1]))
    Y_M.append(float(row[2]))
    Z_M.append(float(row[3]))
f.close()

sma = 383717

RE = 6371.1
RM = 1731.1

# plt.plot(T)

fig = plt.figure()
ax3d = plt.axes(projection='3d')

ax3d.scatter(X[0],Y[0],Z[0],c='m')
ax3d.plot(X,Y,Z,'m')
ax3d.plot(X_E,Y_E,Z_E,'c')
ax3d.plot(X_M,Y_M,Z_M,'k')
ax3d.scatter(X_M[0],Y_M[0],Z_M[0],c='k')
ax3d.set_xlim([-0.5*sma,0.5*sma])
ax3d.set_ylim([-0.5*sma,0.5*sma])
ax3d.set_zlim([-0.5*sma,0.5*sma])

CS = [ [-0.681734548418712, 0.72788989999491, 0.0735819201824379],[-0.726058751056763, -0.685495600222748, 0.0541707678471442],[0.0898704373319078, -0.0164947131122515, 0.995816865158005]]
origin = [-1.32523475880271, 1.41495688941222, 0.143037078685438]
color = ['r','g','b']
scale = 50000
for i in range(3):
    ax3d.plot([origin[0], origin[0] + scale*CS[i][0]],[origin[1], origin[1] + scale*CS[i][1]],[origin[2], origin[2] + scale*CS[i][2]],color[i])


fig = plt.figure()
ax3d = plt.axes(projection='3d')

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xS = np.cos(u)*np.sin(v)
yS = np.sin(u)*np.sin(v)
zS = np.cos(v)

xE = RE*xS - mu*sma
xM = RM*xS + mu1*sma
ax3d.plot_wireframe(xE, RE*yS, RE*zS, color="b")
# ax3d.plot_wireframe(xM, RM*yS, RM*zS, color="r")
ax3d.scatter(L1*sma,0,0,marker='*',color="y")

ax3d.plot(X_EML1,Y_EML1,Z_EML1,'r')
ax3d.set_xlim([(L1 - 0.2)*sma,(L1 + 0.2)*sma])
ax3d.set_ylim([-0.2*sma,0.2*sma])
ax3d.set_zlim([-0.2*sma,0.2*sma])

plt.show()

