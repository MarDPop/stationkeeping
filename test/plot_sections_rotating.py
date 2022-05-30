import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import dynamics

mu = 0.01215058426954
mu1 = 1-mu
mean_motion = 2.6590930417337446e-06
L1 = dynamics.getL1(mu)

X = []
Y = []
Z = []
T = []
for i in range(30):
    f = open("../output/section_rotating_" + str(i), "r")
    data = f.readlines()

    for line in data:
        row = line.split()
        T.append(float(row[0]))
        X.append(float(row[1]))
        Y.append(float(row[2]))
        Z.append(float(row[3]))

    f.close()


sma = 383717

RE = 6371.1
RM = 1731.1

fig = plt.figure()
ax3d = plt.axes(projection='3d')

ax3d.plot(X,Y,Z,'m')
ax3d.set_xlim([-0.5*sma,0.5*sma])
ax3d.set_ylim([-0.5*sma,0.5*sma])
ax3d.set_zlim([-0.5*sma,0.5*sma])

plt.show()

exit(1)

CS = [ [1, 0, 0],[0, 1, 0],[0, 0, 1]]
origin = [0,0,0]
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
ax3d.plot_wireframe(xM, RM*yS, RM*zS, color="r")
ax3d.scatter(L1*sma,0,0,marker='*',color="y")

plt.show()

