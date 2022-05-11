import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import dynamics

mu = 0.01215058426954
mu1 = 1-mu
mean_motion = 2.6590930417337446e-06
sma = 385000
L1 = dynamics.getL1(mu)

f = open("test_orbit", "r")
data = f.readlines()

X = []
Y = []
Z = []
T = []
for line in data:
    row = line.split()
    T.append(float(row[0]))
    X.append(float(row[1])*sma)
    Y.append(float(row[2])*sma)
    Z.append(float(row[3])*sma)
f.close()

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

