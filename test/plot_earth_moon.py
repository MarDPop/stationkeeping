import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import dynamics

mu = 0.01215058426954
mu1 = 1-mu
mean_motion = 2.6590930417337446e-06
L1 = dynamics.getL1(mu)

f = open("../output/earth_orbit", "r")
data = f.readlines()

Xe = []
Ye = []
Ze = []
Te = []
for line in data:
    row = line.split()
    Te.append(float(row[0]))
    Xe.append(float(row[1]))
    Ye.append(float(row[2]))
    Ze.append(float(row[3]))

f.close()

f = open("../output/moon_orbit", "r")
data = f.readlines()

Xm = []
Ym = []
Zm = []
Tm = []
for line in data:
    row = line.split()
    Tm.append(float(row[0]))
    Xm.append(float(row[1]))
    Ym.append(float(row[2]))
    Zm.append(float(row[3]))
    
f.close()

fig = plt.figure()
ax3d = plt.axes(projection='3d')

ax3d.plot(Xe,Ye,Ze,'b')
ax3d.plot(Xm,Ym,Zm,'k')

plt.show()

