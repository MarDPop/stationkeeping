# -*- coding: utf-8 -*-
from tkinter import *
# from tkinter.ttk import *
import computations
import plot
import dynamics
import solarsystem
import subprocess
import numpy as np

class App:
    def __init__(self, win):
        
        col1 = 8
        col2 = 216
        rowHeight = 32
        
        row = 8
        lbl = Label(win, text='Circular Restricted Halo Orbit Parameters',font='Helvetica 12 bold')
        lbl.place(x=col1, y=row)
        
        self.computation = IntVar()
        self.computation.set(1)
        row += rowHeight
        self.rIdeal = Radiobutton(win, text="Circular Restricted", variable=self.computation,value=1)
        self.rIdeal.place(x=10,y=row)
        self.rMotions = Radiobutton(win, text="Lunar Motion", variable=self.computation,value=3)
        self.rMotions.place(x=180,y=row)
        
        row += rowHeight+10
        lbl = Label(win, text='Start:')
        lbl.place(x=8, y=row)
        self.hour = Entry(win,width=3)
        self.hour.place(x=60,y=row)
        self.hour.insert(0,"HH")
        
        lbl = Label(win, text=':')
        lbl.place(x=95, y=row)
        
        self.min = Entry(win,width=3)
        self.min.place(x=110,y=row)
        self.min.insert(0,"MM")
        
        lbl = Label(win, text=':')
        lbl.place(x=145, y=row)
        
        self.sec = Entry(win,width=3)
        self.sec.place(x=160,y=row)
        self.sec.insert(0,"SS")
        
        lbl = Label(win, text='(UTC)')
        lbl.place(x=195, y=row)
        
        self.hour = Entry(win,width=3)
        self.hour.place(x=250,y=row)
        self.hour.insert(0,"DD")
        
        lbl = Label(win, text='/')
        lbl.place(x=285, y=row)
        
        self.min = Entry(win,width=3)
        self.min.place(x=300,y=row)
        self.min.insert(0,"MM")
        
        lbl = Label(win, text='/')
        lbl.place(x=335, y=row)
        
        self.sec = Entry(win,width=5)
        self.sec.place(x=350,y=row)
        self.sec.insert(0,"YYYY")
        
        row += rowHeight
        self.solarPressure = Checkbutton(win, text="Add Sun")
        self.solarPressure.place(x=col1,y=row)
        self.solarPressure.deselect()
        
        row += rowHeight
        self.solarPressure = Checkbutton(win, text="Add Solar Radiation Pressure")
        self.solarPressure.place(x=col1,y=row)
        self.solarPressure.deselect()
        
        row += rowHeight
        lbl = Label(win, text='Satellite Cross Area (m2):')
        lbl.place(x=col1, y=row)
        self.satArea = Entry()
        self.satArea.place(x=col2,y=row)
        self.satArea.insert(0,"0.5")
        
        row += rowHeight
        lbl = Label(win, text='Satellite Mass (kg):')
        lbl.place(x=col1, y=row)
        self.satArea = Entry()
        self.satArea.place(x=col2,y=row)
        self.satArea.insert(0,"200")
        
        row += rowHeight + 10
        lbl = Label(win, text='Amplitude Z (km): ')
        lbl.place(x=col1, y=row)
        self.amplitude = Entry()
        self.amplitude.place(x=col2,y=row)
        self.amplitude.insert(0,"10000")
        
        row += rowHeight
        lbl = Label(win, text='Tolerance (m/s):')
        lbl.place(x=col1, y=row)
        self.tolerance = Entry()
        self.tolerance.place(x=col2,y=row)
        self.tolerance.insert(0,"1e-4")
        
        row += rowHeight
        lbl = Label(win, text='Earth Moon Distance (km):')
        lbl.place(x=col1, y=row)
        self.distance = Entry()
        self.distance.place(x=col2,y=row)
        self.distance.insert(0,"385000")
        
        row += rowHeight
        lbl = Label(win, text='Integrate Time (days):')
        lbl.place(x=col1, y=row)
        self.days2plot = Entry()
        self.days2plot.place(x=col2,y=row)
        self.days2plot.insert(0,"2")
        
        row += rowHeight
        b1 = Button(win, text='Generate Halo', command=self.computeHalo, bg='#333', fg='White')
        b1.place(x=10, y=row)
        b1 = Button(win, text='Generate Vertical', command=self.computeVertical, bg='#333', fg='White')
        b1.place(x=120, y=row)
        b1 = Button(win, text='Integrate', command=self.plotOrbit, bg='#333', fg='White', width=7)
        b1.place(x=250, y=row)
        b1 = Button(win, text='Generate Heat Map', command=self.createHeatMap, bg='#333', fg='White')
        b1.place(x=320, y=row)

        row += rowHeight+10
        lbl = Label(win, text='Station Keeping Calculation Flags',font='Helvetica 12 bold')
        lbl.place(x=col1, y=row)
        
        row += rowHeight
        lbl = Label(win, text='Evaluation Duration (days):')
        lbl.place(x=8, y=row)
        self.runTime = Entry(win)
        self.runTime.place(x=col2,y=row)
        self.runTime.insert(0,"60")
        
        row += rowHeight
        lbl = Label(win, text='Target Recovery (days):')
        lbl.place(x=8, y=row)
        self.goalTime = Entry(win)
        self.goalTime.place(x=col2,y=row)
        self.goalTime.insert(0,"15")
        
        row += rowHeight
        self.usePertubation = IntVar()
        perturbations = Checkbutton(win, text="Add Random Perturbation", variable=self.usePertubation)
        perturbations.place(x=col1,y=row) 
        
        row += rowHeight
        lbl = Label(win, text='Perturbation 1-sigma (mm/s2):')
        lbl.place(x=col1, y=row)
        self.perturbSize = Entry()
        self.perturbSize.place(x=col2,y=row)
        self.perturbSize.insert(0,"0.01")

        row += rowHeight
        lbl = Label(win, text='Perturbation Interval (days):')
        lbl.place(x=col1, y=row)
        self.perturbInterval = Entry()
        self.perturbInterval.place(x=col2,y=row)
        self.perturbInterval.insert(0,"5")
  
        row += rowHeight
        b1 = Button(win, text='Compute Yearly Delta-V', command=self.computeDV, bg='#333', fg='White')
        b1.place(x=10, y=row)
        
        row += rowHeight+10
        self.yearlyDv = Label(win, text='? m/s',font='Helvetica 12 bold')
        self.yearlyDv.place(x=col1, y=row)
        
        self.initial_state = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        self.orbit_times = []
        self.orbit_X = []
        self.orbit_Y = []
        self.orbit_Z = []
        
    def computeHalo(self):
        Az0 = float(self.amplitude.get())
        tol = float(self.tolerance.get())*1e-3
        sma = float(self.distance.get())
        
        cmd = "./DC.exe -t " + str(tol) + " -e " + str(Az0) + " -r " + str(sma)
        
        print(cmd)
        
        proc = subprocess.Popen(cmd, shell=False, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out = proc.communicate()[0]
        
        f = open("output", "r")
        data = f.read()
        print("initial state: " + data)
        arr = data.split()
        print(len(arr))
        for idx in range(6):
            self.initial_state[idx] = float(arr[idx])

        self.computeOrbit()

    def computeOrbit(self):
        print('computing orbit')
        mu = 0.012155650403206974
        mu1 = 0.987844349596793
        mean_motion = 2.6590930417337446e-06
        
        x = self.initial_state.copy()
        
        t = 0
        t_record = 0
        record_interval = 0.01
        dt = 1e-5
        t_final = 100

        # RK Constants
        dt2 = dt/2
        dt3 = dt/3
        dt6 = dt/6
        dt8 = dt/8

        self.orbit_times = []
        self.orbit_X = []
        self.orbit_Y = []
        self.orbit_Z = []
        while(t < t_final):
            xn = x.copy()
            
            if t > t_record:
                self.orbit_X.append(xn[0])
                self.orbit_Y.append(xn[1])
                self.orbit_Z.append(xn[2])
                self.orbit_times.append(t)
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
            
            if x[1] > 0 and xn[1] < 0 and t > 0.01:
                dydt = (x[1] - xn[1])/dt
                dt2 = -xn[1]/dydt

                dxdt = (x[0] - xn[0])/dt
                dzdt = (x[2] - xn[2])/dt
                self.orbit_X.append(xn[0] + dxdt*dt2)
                self.orbit_Y.append(xn[1] + dydt*dt2)
                self.orbit_Z.append(xn[2] + dzdt*dt2)
                self.orbit_times.append(t + dt2)
                break

            t += dt

        days = self.orbit_times[-1]/mean_motion/86400
        print("Orbit period: " + "%5.2f" % days)
        
    def computeVertical(self):
        Az0 = float(self.amplitude.get())
        nIter = int(self.nIterations.get())
        tol = float(self.tolerance.get())*1e-3
        sma = float(self.distance.get())
        
        x0 = computations.initialStateEstimate(Az0,sma)
        
        self.initial_state = computations.DC(x0, nIter, tol,sma)
        self.position.delete(0,END)
        p = np.array2string(self.initial_state[:3], precision=12, separator=',',suppress_small=True)
        self.position.insert(0,p)
        
        self.velocity.delete(0,END)
        v = np.array2string(self.initial_state[3:], precision=12, separator=',',suppress_small=True)
        self.velocity.insert(0,v)
        
    def plotOrbit(self):
        t = float(self.days2plot.get())*86400
        sma = float(self.distance.get())
        plot.plots(self.initial_state,t,sma)

    def createHeatMap(self):
        sma = float(self.distance.get())
        plot.heatmap(self.orbit_times, self.orbit_X, self.orbit_Y, self.orbit_Z,sma)
        
    def computeDV(self):
        Az0 = float(self.amplitude.get())
        sma = float(self.distance.get())
        runtime = float(self.runTime.get())*86400
        t_t = float(self.goalTime.get())*86400
        if self.usePertubation.get():
            perturbation = float(self.perturbSize.get())*1e-6
            perturbTime = float(self.perturbSize.get())*86400
        else:
            perturbation = 0
            perturbTime = -1
        
        ans = computations.computeDeltaV(self.initial_state,t_t,runtime,sma,perturbation,perturbTime)
        
        yearFactor = 365.0/float(self.runTime.get())
        self.yearlyDv.config( text = '%5.1f' % (ans[0]*yearFactor) + " m/s")
        
        plot.plotstationkeeping(ans[3], ans[4], ans[5], ans[1], ans[2], sma)
        
        


window = Tk()
app = App(window)
window.title("Lagrange Orbit Gui")
window.geometry("480x800+10+10")
window.mainloop()
