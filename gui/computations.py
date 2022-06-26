# -*- coding: utf-8 -*-
import numpy as np
import dynamics
import solarsystem

def DC(x0,max_iterations,tol,sma):
    
    mu = solarsystem.MOON_MU/(solarsystem.MOON_MU+solarsystem.EARTH_MU)
    mu1 = 1-mu
    mean_motion = np.sqrt((solarsystem.MOON_MU + solarsystem.EARTH_MU)/(sma**3))
    f = sma*mean_motion
    tol /= f
    
    # dt
    dt = 2e-5
    dtSTM = 0.001
    # RK Constants
    dt2 = dt/2
    dt3 = dt/3
    dt6 = dt/6
    dt8 = dt/8
    
    # final time to integrate to period
    t_final = 5
    xs = []
    ts = []
    # Other variables
    t_record = 0
    it = 0
    alternate = True
    while it < max_iterations:
        t = 0
        x = x0.copy()
        # Runge Kutta Propagate
        xs = []
        ts = []
        t_record = 0
        while(t < t_final):
            xn = x.copy()
            
            if t > t_record:
                xs.append(xn)
                ts.append(t)
                t_record += dtSTM
                
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
    
            # Ensure Boundary Condition
            if x[1] < 0:
                delta = -xn[1]/(x[1] - xn[1])
                
                for j in range(6):
                    x[j] = xn[j] + delta*(x[j] - xn[j])
                
                xs.append(x)
                ts.append(t + dt*delta)
                break
    
            t += dt
    
        # Stop Condition
        print("half period x velocity (m/s): " + str(x[3]*f*1e3))
        print("half period z velocity (m/s): " + str(x[5]*f*1e3))
    
        e = (abs(x[3]) + abs(x[5]))
        
        if e < tol:
            break
    
        # Integrate State Transition Matrix
        STM = np.eye(6)
        n = len(xs)-1
        for i in range(n):
            A = dynamics.getA(xs[i],mu,mu1)
            dSTM = A.dot(STM)
            STM = STM + dSTM*(ts[i+1]-ts[i])
    
        # Solve for delta intial state
        dx = dynamics.dxdt(x,t,mu,mu1)
        
        xdd = dx[3]/x[4]
        zdd = dx[5]/x[4]
    
        stepSize = 1/(1+e*5) # step limiter
    
        # alternate between dx0 and dx0 matrices
        if alternate:
            A00 = STM[3][0] - xdd*STM[1][0]
            A01 = STM[3][4] - xdd*STM[1][4]
            A10 = STM[5][0] - zdd*STM[1][0]
            A11 = STM[5][4] - zdd*STM[1][4]
            det = A00*A11 - A10*A01
            dx0 = (A01*x[5] - A11*x[3])/det
            dyd0 = (A10*x[3] - A00*x[5])/det
            x0[0] += stepSize*dx0
        else:
            A00 = STM[3][2] - xdd*STM[1][2]
            A01 = STM[3][4] - xdd*STM[1][4]
            A10 = STM[5][2] - zdd*STM[1][2]
            A11 = STM[5][4] - zdd*STM[1][4]
            det = A00*A11 - A10*A01
            dz0 = (A01*x[5] - A11*x[3])/det
            dyd0 = (A10*x[3] - A00*x[5])/det
            x0[2] += stepSize*dz0
    
        x0[4] += stepSize*dyd0
    
        print("corrected initial state")
        print(x0)
    
        alternate = not alternate
        
        it += 1
        
    return x0

def initialStateEstimate(Az0,sma):

    mean_motion = np.sqrt((solarsystem.EARTH_MU + solarsystem.MOON_MU)/(sma**3))
    mu = solarsystem.MOON_MU/(solarsystem.MOON_MU + solarsystem.EARTH_MU)
    
    # Get Initial State Estimate
    L1 = dynamics.getL1(mu)
    
    x0 = dynamics.getInitialState(L1,mu,mean_motion,sma,Az0)
    print("initial 3rd order state guess")
    print(x0)
    
    return x0

def rkstep(x,t,dt):
    dt3 = dt*0.33333333333333333
    dt6 = dt*0.16666666666666666
    dt2 = dt*0.5
    xn = x.copy()
    k0 = dynamics.dxdt(x,t,0.012155650403206974,0.987844349596793) 
    x = xn + k0*dt3
    k1 = dynamics.dxdt(x,t + dt3,0.012155650403206974,0.987844349596793)
    x = xn + (k0 + k1)*dt6
    k2 = dynamics.dxdt(x,t + dt3,0.012155650403206974,0.987844349596793)*3
    x = xn + (k0 + k2)*(dt*0.125)
    k3 = dynamics.dxdt(x,t + dt2,0.012155650403206974,0.987844349596793)*4
    x = xn + (k0 - k2 + k3)*dt2
    k4 = dynamics.dxdt(x,t + dt,0.012155650403206974,0.987844349596793)
    x = xn + (k0 + k3 + k4)*dt6
    return x

def computeDeltaV(x0,t_goal,t_final,sma,perturbation,perturbTime):
    
    mean_motion = np.sqrt((solarsystem.MOON_MU + solarsystem.EARTH_MU)/(sma**3))
    t_goal *= mean_motion
    t_final *= mean_motion
    perturbTime *= mean_motion

    print("Computing Periodic Orbit Positions and STM")    

    # dt
    dt = 10e-6
    dtSTM = 2e-3

    # Other variables
    t_record = 0
    t = 0
    x = x0.copy()

    nominal_x = []
    nominal_t = []
    As = []
    t_record = 0
    while(t < 100):
        
        if t > t_record:
            nominal_x.append(x.copy())
            nominal_t.append(t)
            As.append(dynamics.getA(x,0.012155650403206974,0.987844349596793))
            t_record += dtSTM
            
        y_old = x[1]
        x = rkstep(x,t,dt)
        t += dt
        
        if x[1] > 0 and x[4] > 0 and y_old < 0:
            break
        
    print("Finished, Perturbing Initial State")   
    
    dY = np.array([1-2*np.random.rand(),10-20*np.random.rand(),1-2*np.random.rand(), 1-2*np.random.rand(),10-20*np.random.rand(),1-2*np.random.rand()])
    dY[:3] = dY[:3]*1e-4/sma
    dY[3:] = dY[3:]*1e-4/(sma*mean_motion)
    
    x = x0 + dY
    
    add_perturb = False
    if perturbation > 1e-10:
        v_perturb = perturbation/(sma*mean_motion*mean_motion)*dt
        perturb = np.array([0,0,0,v_perturb*np.random.randn(),v_perturb*np.random.randn(),v_perturb*np.random.randn()])
        add_perturb = True
    
    mindV = 2e-4/(sma*mean_motion)
    mindV *= mindV
    
    mindX = 10/sma
    mindX *= mindX
    
    nI = int(t_goal/dtSTM)
    n = len(nominal_x) - 1
    t = 0
    t_period = 0
    t_record = 0
    t_perturb = 0
    t_min = 0.1
    t_last = 0
    dVs = []
    dVtime = []
    dt = 10e-6
    stationkeeping_x = []
    stationkeeping_t = []
    sumdV = 0
    
    print("Computing Station Keeping Burns. Time to run: " + str(t_final))   
    while(t < t_final):            
        
        if t > t_record:
            if add_perturb and t > t_perturb:
                t_perturb += perturbTime
                perturb = np.array([0,0,0,v_perturb*np.random.randn(),v_perturb*np.random.randn(),v_perturb*np.random.randn()])
            
            t_record += 0.02
            stationkeeping_x.append(x.copy())
            stationkeeping_t.append(t)
            print("Time: " + str(t))
            
            if t > t_last: 
            
                t_past_period = t - t_period
                i = int(t_past_period/dtSTM) 
                if i > n:
                    i = n
                    
                nominal = nominal_x[i]
                err = x - nominal
                dx = err[:3]
                if dx.dot(dx) > mindX:
                    STM = np.eye(6)
                    i2 = i
                    count = 0
                    while count < nI:
                        dSTM = As[i2].dot(STM)
                        STM = STM + dSTM*dtSTM
                        if i2 >= n:
                            i2 = -1
                            
                        count += 1
                        i2 += 1
                        
                    A = STM[:3,:3]
                    B = STM[:3,3:]
                    
                    dv = err[3:]
                    
                    dV = -np.linalg.inv(B).dot(A.dot(dx)) - dv
                    
                    dVmag = dV.dot(dV)
                    
                    print(dVmag)
                    
                    if dVmag > mindV:
                        print("Burn at time: " + str(t))
                        dVs.append(dV)
                        dVtime.append(t)
                        t_last = t + t_min
                        x[3] += dV[0]
                        x[4] += dV[1]
                        x[5] += dV[2]
                        sumdV += np.sqrt(dVmag)
        
        y_old = x[1]
        x = rkstep(x,t,dt)
        
        if add_perturb:
            x = x + perturb
        
        if x[1] > 0 and x[4] > 0 and y_old < 0:
            t_period = t - x[1]/x[4]
            print("Orbit completed: " + str(t_period))

        t += dt
        
    sumdV *= (sma*mean_motion*1e3)
    
    return [sumdV,dVs,dVtime,stationkeeping_x,stationkeeping_t,nominal_x,nominal_t]

    

def computeDeltaV_fullOrbital(x0,Az,solarArea,mass,jd0,perturb):
    pass
