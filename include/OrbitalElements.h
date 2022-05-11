#pragma once

#include "Math.h"
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>

#define TWOPI 6.283185307179586476925286766559
#define DEG2RAD 1.7453292519943295769236907685e-2

class OrbitalElements {
public:
	const double mu;
	std::vector< double > jds;
	std::vector< std::array<double,6> > oes;
	std::vector< std::array<double,6> > d_oe;
	std::vector< double > mean_anomaly;
	std::vector< double > mean_motion;
	int idx;

	static const int LIGHTSPEED; 
	static const double EARTH_MU; 
	static const double MOON_MU; 
	static const double SUN_MU; 
	static const double JUPITER_SYSTEM_MU;
	static const double MARS_SYSTEM_MU;
	static const double ICRF_MU;
	static const double EMBARY_MU;
	static const double SOLAR_IRRADIANCE_AU; 
	static const double SOLAR_PRESSURE_AU;
	static const double AU; 
	static const double UT_TDB_DELTAT;
	
	OrbitalElements( const std::vector<double>& JD, const std::vector< std::array<double,6> >& OE, const double& MU) : mu(MU), jds(JD), oes(OE) {
		const int n = jds.size() - 1;
		this->mean_anomaly.reserve(n);
		for(int i = 0; i <= n; i++){
			std::array<double,6>& oe = this->oes[i];
			this->mean_anomaly[i] = OrbitalElements::meanAnomalyFromTrueAnomaly(oe[5],oe[1]);
		}
		
		for(int i = 0; i < n; i++){
			std::array<double,6> d;
			const std::array<double,6>& lo = this->oes[i];
			const std::array<double,6>& hi = this->oes[i+1];
			double djd = this->jds[i+1] - this->jds[i];
			for(i = 0; i < 6;i++) {
				d[i] = (hi[i] - lo[i])/djd;
			}
			this->d_oe.push_back(d);
			this->mean_motion.push_back((this->mean_anomaly[i+1] - this->mean_anomaly[i])/djd);
		}
	}
	
	OrbitalElements(const std::string& fn, const double& MU) : mu(MU){
		
		std::vector< std::string > lines;
		
		std::ifstream myfile;
		myfile.open(fn);

		if(!myfile.is_open()) {
			perror("Error opening file");
			throw 5;
		}
		
		std::string line;
		while(std::getline(myfile, line)) {
			lines.push_back(line);
		}
		myfile.close();
		
		while(lines.back().length() < 10) {
			lines.pop_back();
		}

		lines[0] = lines[0].substr(3); // correct first line;
		
		uint_fast32_t nLines = lines.size();
		uint_fast32_t i = 0;
		uint_fast32_t entry = 0;
		uint_fast32_t entries = (nLines - 4)/5;
		this->jds.resize(entries);
		this->mean_motion.resize(entries);
		this->mean_anomaly.resize(entries);
		this->oes.resize(entries);
		while(entry < entries) {
			this->jds[entry] = std::stod(lines[i].substr(0,18));
			std::array<double,6>& oe = this->oes[entry];
			i++;
			oe[1] = std::stod(lines[i].substr(4,27));
			oe[2] = std::stod(lines[i].substr(57))*DEG2RAD;
			i++;
			oe[3] = std::stod(lines[i].substr(4,27))*DEG2RAD;
			oe[4] = std::stod(lines[i].substr(31,53))*DEG2RAD;
			i++;
			this->mean_motion[entry] = std::stod(lines[i].substr(4,27))*(DEG2RAD*86400); // for some reason HORIZONS mean motion is off
			this->mean_anomaly[entry] = std::stod(lines[i].substr(31,53))*DEG2RAD;
			oe[5] = std::stod(lines[i].substr(57))*DEG2RAD;
			i++;
			oe[0] = std::stod(lines[i].substr(4,27));
			i++;

			entry++;
		}
		entries--;
		for(entry = 0; entry < entries; entry++){
			std::array<double,6> d;
			uint_fast32_t next = entry+1;
			const std::array<double,6>& lo = this->oes[entry];
			const std::array<double,6>& hi = this->oes[next];
			
			double djd = this->jds[next] - this->jds[entry];
			for(i = 0; i < 6;i++) {
				d[i] = (hi[i] - lo[i])/djd;
			}
			if(d[5] < 0) {
				d[5] += TWOPI;
			}
			
			this->d_oe.push_back(d);
			double dMA = this->mean_anomaly[next] - this->mean_anomaly[entry];
			if(dMA < 0){
				dMA += TWOPI;
			}
			this->mean_motion[entry] = dMA/djd;
		}
		
	}
	~OrbitalElements(){}
	
	double getMu() const{
		return this->mu;
	}
	
	static double trueAnomalyFromEccentricAnomaly( const double& EA, const double& eccentricity ){
		double beta = eccentricity/(1 + sqrt(1 - eccentricity*eccentricity));
		double s = sin(EA);
        return EA + 2*atan(beta*s/(1 - beta*sqrt(1 - s*s)));
    }

    static double trueAnomalyFromMeanAnomaly( const double& MA, const double& eccentricity ){
        return trueAnomalyFromEccentricAnomaly(eccentricAnomalyFromMeanAnomaly(MA,eccentricity),eccentricity);
    }
	
	static double trueAnomalyFromMeanAnomalyApprox( const double& MA, const double& eccentricity ){
		double e2 = eccentricity*eccentricity;
        return MA + eccentricity*(2 - 0.25*e2)*sin(MA) + e2*(1.2*sin(MA+MA) + eccentricity*1.083333333333333*sin(3*MA));
    }
    
    static double meanAnomalyFromTrueAnomalyApprox( const double& f, const double& eccentricity){
        double e2 = eccentricity*eccentricity;
        double e3 = e2*eccentricity;
        double f2 = f+f;
        return f + eccentricity*(-2*sin(f) + (0.75*eccentricity+0.125*e3)*sin(f2) - 0.333333333333333*e2*sin(f2+f)+ 0.15625*e3*sin(f2+f2));
    }

    static double eccentricAnomalyFromMeanAnomaly( const double& MA, const double& eccentricity){
        double EA = MA;
        for (int iter = 0; iter < 10; iter++) {
            double dEA = (EA - eccentricity*sin(EA) - MA)/(1 - eccentricity*cos(EA));
            EA -= dEA;
            if (fabs(dEA) < 1e-10)
                break;
        }
        return EA;
		/*
		double f = MA/3.1415926535897932384626433832795;
		int rotations = (int)f;
		bool flip = (rotations % 2 == 1);
		double Ei;
		if(flip){
			Ei = 3.1415926535897932384626433832795*(1 - (f - rotations));
		} else {
			Ei = (f - rotations)*3.1415926535897932384626433832795;
		}
		double e2 = eccentricity*eccentricity;
		for (int iter = 0; iter < 20; iter++) {
			double c = eccentricity*cos(Ei);
			double dE = (Ei - sqrt(e2 - c*c) - MA)/(1 - c);
			Ei -= dE;
			if (fabs(dE) < 1e-10)
				break;
		}
		if(flip) {
			return -Ei;
		} 
		return Ei;
		*/
    }
    
    static double meanAnomalyFromEccentricAnomaly( const double& EA, const double& eccentricity) {
        return EA - eccentricity*sin(EA);
    }
    
    static double eccentricAnomalyFromTrueAnomaly( const double& TA, const double& eccentricity ) {
        return atan2(sqrt(1 - eccentricity*eccentricity)*sin(TA), eccentricity + cos(TA));
    }
    
    static double meanAnomalyFromTrueAnomaly( const double& f, const double& eccentricity){
        return meanAnomalyFromEccentricAnomaly(eccentricAnomalyFromTrueAnomaly(f,eccentricity),eccentricity);
    }
    
    static std::array<double,6> kepler2cartesian(const std::array<double,6>& oe, const double& mu){   
		   
        double tmp = 1.0 - oe[1]*oe[1];
        double st = sin(oe[5]);
        double ct = cos(oe[5]);
        double tmp2 = 1.0 + oe[1]*ct;
        double radius = oe[0]*tmp/tmp2;
        double x = radius*ct;
        double y = radius*st;
        tmp = sqrt( mu*oe[0]*tmp )/(radius*tmp2);
        double v_x = -st*tmp;
        double v_y = (oe[1]+ct)*tmp;
            
        if (fabs(oe[2]) < 1e-8){
			std::array<double,6> out = {x,y,0,v_x,v_y,0};
            return out;
        }
            
        double cw = cos(oe[4]);
        double sw = sin(oe[4]);
        double co = cos(oe[3]);
        double so = sin(oe[3]);
        
        st = sin(oe[2]);
        ct = sqrt(1.0-st*st);
        double Rxx = cw*co - sw*ct*so;
        double Rxy = -(sw*co + cw*ct*so);
        double Ryx = cw*so + sw*ct*co;
        double Ryy = cw*ct*co - sw*so;
        double Rzx = sw*st;
        double Rzy = cw*st;
        
        std::array<double,6> out;
        out[0] = Rxx*x + Rxy*y;
        out[1] = Ryx*x + Ryy*y;
        out[2] = Rzx*x + Rzy*y;
        out[3] = Rxx*v_x + Rxy*v_y;
        out[4] = Ryx*v_x + Ryy*v_y;
        out[5] = Rzx*v_x + Rzy*v_y;
        return out;
    }
	
	static std::array<double,3> kepler2position(const std::array<double,6>& oe){   
        double st = sin(oe[5]);
        double ct = cos(oe[5]);
        double radius = oe[0]*(1.0 - oe[1]*oe[1])/(1.0 + oe[1]*ct);
        double x = radius*ct;
        double y = radius*st;
            
        if (fabs(oe[2]) < 1e-8){
			std::array<double,3> out = {x,y,0};
            return out;
        }
            
        double cw = cos(oe[4]);
        double sw = sin(oe[4]);
        double co = cos(oe[3]);
        double so = sin(oe[3]);
        
        st = sin(oe[2]);
        ct = sqrt(1.0-st*st);
        double Rxx = cw*co - sw*ct*so;
        double Rxy = -(sw*co + cw*ct*so);
        double Ryx = cw*so + sw*ct*co;
        double Ryy = cw*ct*co - sw*so;
        double Rzx = sw*st;
        double Rzy = cw*st;
        
        std::array<double,3> out;
        out[0] = Rxx*x + Rxy*y;
        out[1] = Ryx*x + Ryy*y;
        out[2] = Rzx*x + Rzy*y;
        return out;
    }
    
    static std::array<double,6> cartesian2kepler(const std::array<double,6>& state, const double& mu) {
		std::array<double,6> oe;
        std::array<double,3> h = {state[1]*state[5] - state[2]*state[4], state[2]*state[3] - state[0]*state[5], state[0]*state[4] - state[1]*state[3]};
        
        std::array<double,2> n = {h[1],-h[0]}; // z is implicit 0
        double v2 = state[3]*state[3] + state[4]*state[4] + state[5]*state[5];
        double r_inv = 1/sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
        double rv = state[0]*state[3] + state[1]*state[4] + state[2]*state[5];
        std::array<double,3> e;
        double tmp1 = v2/mu - r_inv;
        double tmp2 = rv/mu;
        e[0] = state[0]*tmp1 + state[3]*tmp2;
        e[1] = state[1]*tmp1 + state[4]*tmp2;
        e[2] = state[2]*tmp1 + state[5]*tmp2;
        double egy = v2/2 - mu*r_inv;
    
        oe[0] = -mu/(2*egy);
        oe[1] = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
        double nmag = n[0]*n[0]+n[1]*n[1];
        oe[2] = acos(h[2]/sqrt(nmag + h[2]*h[2]));
        
        if ( fabs(oe[2]) > 1e-9){
            nmag = sqrt(nmag);
            oe[3] = acos(n[0]/nmag);
            if (n[1] < 0){
                oe[3] = TWOPI - oe[3];
            }
            	
            oe[4] = acos((n[0]*e[0]+n[1]*e[1])/(nmag*oe[1]));
            if (e[2] < 0){
                oe[4] = TWOPI - oe[4];
            }
        }
        
        oe[5] = acos((e[0]*state[0]+e[1]*state[1]+e[2]*state[2])*r_inv/oe[1]);
        if (rv < 0){
            oe[5] = TWOPI - oe[5];
        }
        return oe;
    }
	
	std::array<double,6> get(const double& jd) const {
		
		if(jd < this->jds[0]){
			return this->oes[0];
		}

		if(jd > this->jds.back()){
			return this->oes.back();
		}
		
		std::array<double,6> oe_interp;
		uint_fast32_t idx = Math::lowerbound(this->jds.data(),jd,this->jds.size());
		
		double djd = jd - this->jds[idx];
		const std::array<double,6>& oe = this->oes[idx];
		const std::array<double,6>& slope = this->d_oe[idx];
		for(int i = 0; i < 6; i++) {
			oe_interp[i] = oe[i] + slope[i]*djd;
		}
		return oe_interp;
	}
	
	std::array<double,6> getSimple(const double& jd) {
		if(jd > this->jds.back()){
			return this->oes.back();
		}
		
		while(jd > this->jds[this->idx+1]){
			this->idx++;
		}
		
		while(jd < this->jds[this->idx] && this->idx > 0){
			this->idx--;
		}
		
		std::array<double,6> oe_interp;
		double djd = jd - this->jds[this->idx];
		const std::array<double,6>& oe = this->oes[this->idx];
		const std::array<double,6>& slope = this->d_oe[this->idx];
		for(int i = 0; i < 6; i++) {
			oe_interp[i] = oe[i] + slope[i]*djd;
		}
		return oe_interp;
	}
	
};

 const int OrbitalElements::LIGHTSPEED = 299792458; // m/s
 // updated values from JPL
 // https://ssd.jpl.nasa.gov/astro_par.html
 const double OrbitalElements::EARTH_MU = 398600.435507; // km3 / s2 
 const double OrbitalElements::MOON_MU = 4902.800118; // km3 / s2
 const double OrbitalElements::SUN_MU =  1.32712440018e11; // km3 / s2 
 const double OrbitalElements::JUPITER_SYSTEM_MU = 126712764.100000; 
 const double OrbitalElements::MARS_SYSTEM_MU = 42828.375816;
 const double OrbitalElements::ICRF_MU = 1.3288932302018904e11;
 const double OrbitalElements::EMBARY_MU = 4.0350323562548013e05;
 
 const double OrbitalElements::SOLAR_IRRADIANCE_AU = 1380; // W/m2 at 1 AU
 const double OrbitalElements::SOLAR_PRESSURE_AU = 4.6031845137344982841429586597539e-6;
 const double OrbitalElements::AU = 149597870.7; // km
 const double OrbitalElements::UT_TDB_DELTAT = 69.185602; // value at 4/20/22 will be +/- 0.005 sec for next decade