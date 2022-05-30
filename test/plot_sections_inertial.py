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
    f = open("../output/section_inertial_" + str(i), "r")
    data = f.readlines()

    for line in data:
        row = line.split()
        T.append(float(row[0]))
        X.append(float(row[1]))
        Y.append(float(row[2]))
        Z.append(float(row[3]))

    f.close()

f = open("../output/earth_orbit", "r")
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

f = open("../output/moon_orbit", "r")
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

fig = plt.figure()
ax3d = plt.axes(projection='3d')

ax3d.plot(X,Y,Z,'m')
ax3d.plot(X_E,Y_E,Z_E,'c')
ax3d.plot(X_M,Y_M,Z_M,'k')
ax3d.set_xlim([-0.5*sma,0.5*sma])
ax3d.set_ylim([-0.5*sma,0.5*sma])
ax3d.set_zlim([-0.5*sma,0.5*sma])

plt.show()
