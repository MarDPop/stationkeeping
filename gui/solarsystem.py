# -*- coding: utf-8 -*-
import numpy as np

# CONSTANTS
LIGHTSPEED =  299792458 # m/s
EARTH_MU = 3.986004418e5 # km2 / s3
MOON_MU = 4.9048695e3 # km2 / s3
SUN_MU = 1.179e11 # km2 / s3
ICRF_MU = 1.3288932302018904E+11
EMBARY_MU = 4.0350323562548013E+05
SOLAR_IRRADIANCE_AU = 1380 # W/m2 at 1 AU
SOLAR_PRESSURE_AU = SOLAR_IRRADIANCE_AU/LIGHTSPEED
AU = 149597870.7 # km

sunJD = []
sunICRF = []
earthJD = []
earthICRF = []
moonJD = []
moonECI = []

def loadData():
    sunFile = open('sun_ICRF.txt','r')
    
    lines = sunFile.readlines()

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if 'TDB' in line:
            jd_str = line.split('=')
            jd = float(jd_str[0].strip()[3:])
            i += 1
            line = lines[i]
            ecc = float(line[4:27])
            inc = float(line[57:])*0.017453292519943295
            i += 1
            line = lines[i]
            omega = float(line[4:27])*0.017453292519943295
            argP = float(line[31:53])*0.017453292519943295
            i += 1
            line = lines[i]
            TA = float(line[57:])*0.017453292519943295
            i += 1
            line = lines[i]
            sma = float(line[4:27])
            oe = np.array([sma,inc,ecc,omega,argP,TA])
            
            sunJD.append(jd)
            sunICRF.append(oe)
            
        i+=1
        
        earthFile = open('earth_ICRF.txt','r')
    
        lines = earthFile.readlines()
    
        i = 0
        n = len(lines)
        while i < n:
            line = lines[i]
            if 'TDB' in line:
                jd_str = line.split('=')
                jd = float(jd_str[0].strip()[3:])
                i += 1
                line = lines[i]
                ecc = float(line[4:27])
                inc = float(line[57:])*0.017453292519943295
                i += 1
                line = lines[i]
                omega = float(line[4:27])*0.017453292519943295
                argP = float(line[31:53])*0.017453292519943295
                i += 1
                line = lines[i]
                TA = float(line[57:])*0.017453292519943295
                i += 1
                line = lines[i]
                sma = float(line[4:27])
                oe = np.array([sma,inc,ecc,omega,argP,TA])
                
                earthJD.append(jd)
                earthICRF.append(oe)
                
            i+=1
        
        moonFile = open('moon_ECI.txt','r')
    
        lines = moonFile.readlines()
    
        i = 0
        n = len(lines)
        while i < n:
            line = lines[i]
            if 'TDB' in line:
                jd_str = line.split('=')
                jd = float(jd_str[0].strip()[3:])
                i += 1
                line = lines[i]
                ecc = float(line[4:27])
                inc = float(line[57:])*0.017453292519943295
                i += 1
                line = lines[i]
                omega = float(line[4:27])*0.017453292519943295
                argP = float(line[31:53])*0.017453292519943295
                i += 1
                line = lines[i]
                TA = float(line[57:])*0.017453292519943295
                i += 1
                line = lines[i]
                sma = float(line[4:27])
                oe = np.array([sma,inc,ecc,omega,argP,TA])
                
                moonJD.append(jd)
                moonECI.append(oe)
                
            i+=1
        
def orbitalElements2cartesian(oe,mu):
    tmp = 1.0 - oe[1]*oe[1]
    st = np.sin(oe[5])
    ct = np.cos(oe[5])
    tmp2 = 1.0 + oe[1]*ct
    radius = oe[0]*tmp/tmp2
    x = radius*ct
    y = radius*st
    tmp = np.sqrt( mu*oe[0]*tmp )/(radius*tmp2)
    v_x = -st*tmp
    v_y = (oe[1]+ct)*tmp
        
    if abs(oe[2]) < 1e-8:
        return np.array([x,y,0,v_x,v_y,0])
        
    cw = np.cos(oe[4])
    sw = np.sin(oe[4])
    co = np.cos(oe[3])
    so = np.sin(oe[3])
    
    st = np.sin(oe[2])
    ct = np.sqrt(1.0-st*st)
    Rxx = cw*co - sw*ct*so
    Rxy = -(sw*co + cw*ct*so)
    Ryx = cw*so + sw*ct*co
    Ryy = cw*ct*co - sw*so
    Rzx = sw*st
    Rzy = cw*st
    
    out = np.zeros((1,6))
    out[0] = Rxx*x + Rxy*y
    out[1] = Ryx*x + Ryy*y
    out[2] = Rzx*x + Rzy*y
    out[3] = Rxx*v_x + Rxy*v_y
    out[4] = Ryx*v_x + Ryy*v_y
    out[5] = Rzx*v_x + Rzy*v_y
    return out

def orbitalElements2position(oe,mu):

    ct = np.cos(oe[5])
    radius = oe[0]*(1.0 - oe[1]*oe[1])/(1.0 + oe[1]*ct)
    x = radius*ct
    y = radius*np.sin(oe[5])
        
    if abs(oe[2]) < 1e-8:
        return np.array([x,y,0])
        
    cw = np.cos(oe[4])
    sw = np.sin(oe[4])
    co = np.cos(oe[3])
    so = np.sin(oe[3])
    
    st = np.sin(oe[2])
    ct = np.sqrt(1.0-st*st)
    Rxx = cw*co - sw*ct*so
    Rxy = -(sw*co + cw*ct*so)
    Ryx = cw*so + sw*ct*co
    Ryy = cw*ct*co - sw*so
    Rzx = sw*st
    Rzy = cw*st
    
    out = np.zeros((3,1))
    out[0] = Rxx*x + Rxy*y
    out[1] = Ryx*x + Ryy*y
    out[2] = Rzx*x + Rzy*y
    return out

def interp(jd,jds,oes):
    if jd <= jds[0]:
        return oes[0]
    
    if jd >= jds[-1]:
        return oes[-1]
            
    lo = 0
    hi = len(jds) - 1
    mid = int((lo+hi)/2)

    while lo != mid:
        if jd > jds[mid]:
            lo = mid
        else:
            hi = mid
            
        mid = int((lo+hi)/2)
        
    oe = oes[lo] + (oes[lo+1] - oes[lo])*(jd - jds[lo])/(jds[lo+1] - jds[lo])
    return oe


def getSunICRF(jd2000):
    oe = interp(jd2000,sunJD,sunICRF)
    return orbitalElements2position(oe,ICRF_MU)

def getEarthICRF(jd2000):
    oe = interp(jd2000,earthJD,earthICRF)
    return orbitalElements2position(oe,ICRF_MU)

def getMoonECI(jd2000):
    oe = interp(jd2000,moonJD,moonECI)
    return orbitalElements2position(oe,EMBARY_MU)

def getSunECI(jd2000):
    earth = getEarthICRF(jd2000)
    sun = getSunICRF(jd2000)
    return sun - earth

def getEarthMoonECICS(jd2000):
    moon_oe = interp(jd2000,moonJD,moonECI)
    state = orbitalElements2cartesian(moon_oe,EMBARY_MU)
    r = state[0:3]
    h = np.cross(r,state[3:])
    z = h / np.sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2])
    x = r / np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
    y = np.cross(z,x)
    
    return np.array([x,y,z])
    