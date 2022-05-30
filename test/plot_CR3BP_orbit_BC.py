import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import dynamics

mu = 0.01215058426954
mu1 = 1-mu
mean_motion = 2.6590930417337446e-06
L1 = dynamics.getL1(mu)

f = open("../output/CR3BP_orbit_converted_BC", "r")
data = f.readlines()

X = []
Y = []
Z = []
T = []
for line in data:
    row = line.split()
    T.append(float(row[0]))
    X.append(float(row[1]))
    Y.append(float(row[2]))
    Z.append(float(row[3]))
    
f.close()

fig = plt.figure()
ax3d = plt.axes(projection='3d')

ax3d.plot(X,Y,Z,'m')
ax3d.scatter(-mu*391000,0,0,'b')
ax3d.scatter(mu1*391000,0,0,'k')

ax3d.set_xlim([-100000,400000])
ax3d.set_ylim([-250000,250000])
ax3d.set_zlim([-250000,250000])

plt.show()

