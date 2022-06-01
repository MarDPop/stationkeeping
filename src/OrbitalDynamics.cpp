#include "../include/OrbitalDynamics.h"

#include "../include/OrbitalElements.h"
#include "../include/Math.h"
#include "../include/Util.h"
#include <cmath>
#include <iostream>
#include <fstream>

std::array<double,3> OrbitalDynamics::convert_cr3bp_to_inertial_pos(EarthMoonSun* dynamics, const std::array<double,6>& state, const double& jd){
	std::array<double,3> r_moon =  dynamics->moon->getPos(jd);
    std::array<double,3> r_earth =  dynamics->earth->getPos(jd);
    std::array<double,3> x = {r_moon[0] - r_earth[0],r_moon[1] - r_earth[1],r_moon[2] - r_earth[2]};
    
    std::array<double,3> v_earth = dynamics->earth->getVel(jd);
    std::array<double,3> h_earth = Math::cross(r_earth,v_earth);
    std::array<double,3> v_moon = dynamics->moon->getVel(jd);
    std::array<double,3> h_moon = Math::cross(r_moon,v_moon);
    
    std::array<double,3> z = {h_earth[0] + h_moon[0],h_earth[1] + h_moon[1],h_earth[2] + h_moon[2]};
    double z_mag = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    double sma = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    for(int i = 0; i < 3; i++) {
        x[i] /= sma;
        z[i] /= z_mag;
    }
    std::array<double,3> y = Math::cross(z,x);

	std::array<double,3> inertial = {state[0]*sma,state[1]*sma,state[2]*sma};

	std::array< std::array<double,3>, 3> CST = {{{x[0],y[0],z[0]},{x[1],y[1],z[1]},{x[2],y[2],z[2]}}};

	return Math::mult(CST,inertial);
}

std::array<double,3> OrbitalDynamics::convert_cr3bp_to_rotating_barycenter(EarthMoonSun* dynamics, const std::array<double,6>& state, const double& jd){
    std::array<double,3> r_moon =  dynamics->moon->getPos(jd);
    std::array<double,3> r_earth =  dynamics->earth->getPos(jd);
    std::array<double,3> x = {r_moon[0] - r_earth[0],r_moon[1] - r_earth[1],r_moon[2] - r_earth[2]};
    double sma = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

	std::array<double,3> inertial = {state[0]*sma,state[1]*sma,state[2]*sma};

	return inertial;
}

std::array<double,6> OrbitalDynamics::convert_cr3bp_to_inertial(EarthMoonSun* dynamics, const std::array<double,6>& state_cr3bp, const double& jd){
	
    std::array<double,3> r_moon = dynamics->moon->getPos(jd);
    std::array<double,3> r_earth = dynamics->earth->getPos(jd);
    std::array<double,3> x = {r_moon[0] - r_earth[0],r_moon[1] - r_earth[1],r_moon[2] - r_earth[2]};
    
    std::array<double,3> v_earth = dynamics->earth->getVel(jd);
    std::array<double,3> h_earth = Math::cross(r_earth,v_earth);
    std::array<double,3> v_moon = dynamics->moon->getVel(jd);
    std::array<double,3> h_moon = Math::cross(r_moon,v_moon);
    
    std::array<double,3> z = {h_earth[0] + h_moon[0],h_earth[1] + h_moon[1],h_earth[2] + h_moon[2]};
    double z_mag = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    double sma = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    for(int i = 0;i < 3; i++) {
        x[i] /= sma;
        z[i] /= z_mag;
    }
    std::array<double,3> y = Math::cross(z,x);

    CR3BP cr3bp(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);
    
    std::array<double,6> state = cr3bp.convert_state_to_inertial(state_cr3bp);

    const double djd = 0.0005;
	
	std::array<double,3> r_moon_2 = dynamics->moon->getPos(jd + djd);
	std::array<double,3> r_earth_2 = dynamics->earth->getPos(jd + djd);
	std::array<double,3> sma_2 = {r_moon_2[0] - r_earth_2[0], r_moon_2[1] - r_earth_2[1],r_moon_2[2] - r_earth_2[2]};
	
	double dt = djd*86400;
	double dL1dt = (Math::norm(sma_2) - sma)*(CR3BP::EARTH_L1 + cr3bp.mu)/dt;
	
	std::array< std::array<double,3>, 3> CST = {{{x[0],y[0],z[0]},{x[1],y[1],z[1]},{x[2],y[2],z[2]}}};
	
	std::array<double,3> pos = {state[0],state[1],state[2]};
	std::array<double,3> vel = {state[3] + dL1dt,state[4],state[5]};
	
	pos = Math::mult(CST,pos);
	vel = Math::mult(CST,vel);

	state[0] = pos[0];
	state[1] = pos[1];
	state[2] = pos[2];
	state[3] = vel[0];
	state[4] = vel[1];
	state[5] = vel[2];
	
	return state;
}

CR3BP::CR3BP() : mu(0.012155650403206974), mu1(0.987844349596793), mu_body1(3.986004418e5) , mu_body2(4.9048695e3) , sma(385000) , mean_motion(2.6590930417337446e-06)  {
}

CR3BP::CR3BP(const double& MU1, const double& MU2, const double& SMA) : mu(MU2/(MU1+MU2)), mu1(MU1/(MU1+MU2)), mu_body1(MU1) , mu_body2(MU2) , sma(SMA), mean_motion(sqrt((MU1 + MU2)/(sma*SMA*SMA))) {
}

void CR3BP::getA(const std::array<double,6>& x, const double& jd, double** A) const{
    Math::clear(A,6);
    A[0][3] = 1;
    A[1][4] = 1;
    A[2][5] = 1;
    A[3][4] = 2;
    A[4][3] = -2;
    
    double dx1 = x[0] + this->mu;
    double dx2 = x[0] - this->mu1;
    double d = x[1]*x[1] + x[2]*x[2];
    double r1_sq = dx1*dx1 + d;
    double r2_sq = dx2*dx2 + d;
    
    double c1 = this->mu1/(r1_sq*sqrt(r1_sq));
    double c2 = this->mu/(r2_sq*sqrt(r2_sq));
    
    double a = c1 + c2;
    
    c1 *= 3/r1_sq;
    c2 *= 3/r2_sq;
    
    double b = c1 + c2;
    
    c1 *= dx1;
    c2 *= dx2;

    d = c1 + c2;
    
    double yz = x[1]*x[2];
    
    A[3][0] = c1*dx1 + c2*dx2 - a + 1;
    A[3][1] = x[1]*d;
    A[3][2] = x[2]*d;
    
    A[4][0] = x[1]*d;
    A[4][1] = x[1]*x[1]*b - a + 1;
    A[4][2] = yz*b;
    
    A[5][0] = x[2]*d;
    A[5][1] = yz*b;
    A[5][2] = x[2]*x[2]*b - a;
}

double CR3BP::getCL1(const double& yL,const int& n){
    double sign = (n % 2 == 0) ? 1 : -1;
    return 1/Math::CUB(yL)* (mu + sign*mu1*std::pow(yL/(1 - yL),n+1));
}

double CR3BP::getCL2(const double& yL,const int& n){
    double sign = (n % 2 == 0) ? 1 : -1;
    return sign/Math::CUB(yL)* (mu + mu1*std::pow(yL/(1 + yL),n+1));
}

double CR3BP::getCL3(const double& yL,const int& n){
    return 1/Math::CUB(yL)* (mu1 + mu*std::pow(yL/(1 + yL),n+1));
}

double CR3BP::getL1(){
    //https://digitalcommons.odu.edu/cgi/viewcontent.cgi?article=1300&context=mae_etds
    //https://map.gsfc.nasa.gov/ContentMedia/lagrange.pdf
    double L1 = 1 - cbrt(mu/3);
    for(int i = 0; i < 10;i++) {
        double r1 = L1 + mu;
        double r2 = L1 - mu1;
        double g1 = mu1/(r1*r1);
        double g2 = mu/(r2*r2);
        double f = L1 - g1 + g2;
        double df = 1 + 2*(g1/r1 - g2/r2);			
        L1 -= 0.7*f/df;
    }
    return L1;
}

double CR3BP::getL2(){
    double L2 = 1 + cbrt(mu/3);
    for(int i = 0; i < 10;i++) {
        double f = L2 - mu1/Math::SQ(L2 + mu) - mu/Math::SQ(L2 - mu1);
        double df = 1 + 2*(mu1/Math::CUB(L2 + mu) + mu/Math::CUB(L2 - mu1));			
        L2 -= 0.7*f/df;
    }
    return L2;
}

double CR3BP::getL3(){
    double L3 = -(1 + 0.4166666*mu);
    for(int i = 0; i < 10;i++) {
        double f = L3 + mu1/Math::SQ(L3 + mu) + mu/Math::SQ(L3 - mu1);
        double df = 1 - 2*(mu1/Math::CUB(L3 + mu) + mu/Math::CUB(L3 - mu1));			
        L3 -= 0.7*f/df;
    }
    return L3;
}

double CR3BP::get_period(const double& Az0) {
    double xL1 = this->getL1();
    double yL = xL1 - this->mu1; // should be negative for L1
    double yLA = fabs(yL);
    double c2 = this->mu/Math::CUB(yLA) + this->mu1/Math::CUB(1 - yLA);
    double c3 = this->getCL1(yLA,3);
    double c4 = this->getCL1(yLA,4);
    double gamma2 = (2 - c2 + sqrt(c2*(9*c2 - 8)))/2;
    double gamma = sqrt(gamma2);
    double k = 2*gamma/(1 + gamma2 - c2);
    double k2 = k*k;
    
    double d1 = 3*gamma2/k*(k*(6*gamma2-1) - 2*gamma);
    
    double a21 = 0.75*c3*(k2 - 2)/(1 + 2*c2);
    double a22 = 0.75*c3/(1 + 2*c2);
    double a23 = -0.75*c3*gamma/(k*d1)*(3*k*k2*gamma - 6*k*(k-gamma) + 4);
    double a24 = -0.75*c3*gamma/(k*d1)*(2 + 3*k*gamma);
    double b21 = -1.5*c3*gamma/d1*(3*k*gamma-4);
    double b22 = 3*c3*gamma/d1;
    double d21 = -c3/(2*gamma2);
    
    double a1 = -1.5*c3*(2*a21 + a23 + 5*d21) - 0.375*c4*(12 - k2);
    double a2 = 1.5*c3*(a24 - 2*a22) + 1.125*c4;
    
    double tmp = 1/(2*gamma*(gamma*(1 + k2) - 2*k));
    double s1 = tmp*(1.5*c3*(2*a21*(k2 - 2) - a23*(k2 + 2) - 2*k*b21) - 0.375*c4*(k2*(3*k2 - 8) + 8) );
    double s2 = tmp*(1.5*c3*(2*a22*(k2 - 2) + a24*(k2 + 2) + 2*k*b22 + 5*d21) + 0.375*c4*(12 - k2));
    
    double l1 = a1 + 2*gamma2*s1;
    double l2 = a2 + 2*gamma2*s2;
    
    double delta = gamma2 - c2;
    
    double r1 = yLA*this->sma;		
    double Az = Az0/r1;
    
    double n1 = this->mean_motion; // sqrt(CUB(yLA)); // normalized, (might be a better way to compute this haha)
    
    double Az2 = Az*Az;
    double Ax2 = -(Az2*l2 + delta)/l1;
    
    double omega2 = s1*Ax2 + s2*Az2;
    double omega = 1 + omega2;
    
    tmp = gamma*omega*n1;
    
    return 6.283185307179586476925286766559/tmp; 
}

std::array<double,6> CR3BP::get_halo_initial_state_3rd_order(const double& Az0, const double& phi, const double& time, const int& n = 1) {
    
    if( n < 0 || n > 4) {
        throw "bad n";
    }
    
    int sign = (n % 2 == 0) ? 1 : -1;
    
    //https://commons.erau.edu/cgi/viewcontent.cgi?article=1565&context=edt
    //https://hal.archives-ouvertes.fr/hal-00312910v2/document
    //https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1980CeMec..22..241R&defaultprint=YES&filetype=.pdf
    //http://www.cds.caltech.edu/archive/help/uploads/wiki/files/39/thurman-worfolk-1996.pdf
    double xL1 = this->getL1();
    double yL = xL1 - this->mu1; // should be negative for L1
    double yLA = fabs(yL);
    double c2 = this->mu/Math::CUB(yLA) + this->mu1/Math::CUB(1 - yLA);
    double c3 = this->getCL1(yLA,3);
    double c4 = this->getCL1(yLA,4);
    double gamma2 = (2 - c2 + sqrt(c2*(9*c2 - 8)))/2;
    double gamma = sqrt(gamma2);
    double k = 2*gamma/(1 + gamma2 - c2);
    double k2 = k*k;
    
    double tmp = -2*c2*c2 + c2 + 1;
    double D1 = 4*gamma2*(4*gamma2 + (c2 - 2)) + tmp;
    double D2 = 9*gamma2*(9*gamma2 + (c2 - 2)) + tmp; 
    double D3 = 2*gamma*(gamma*(1 + k2) - 2*k);
    
    tmp = 0.75*c3/(1 + 2*c2);
    double a21 = tmp*(k2 - 2);
    double a22 = tmp;
    
    tmp = -0.75*gamma*c3/(k*D1);
    double a23 = tmp*(3*k*(k2*gamma - 2*(k - gamma)) + 4);
    double a24 = tmp*(2 + 3*k*gamma);
    
    tmp = 3*c3*gamma/D1;
    double b21 = -0.5*tmp*(3*k*gamma - 4);
    double b22 = tmp;
    
    double d21 = -c3/(2*gamma2);
    tmp = 3/(64*gamma2);
    double d31 = tmp*(4*c3*a24 + c4);
    double d32 = tmp*(4*c3*(a23 - d21) + c4*(4 + k2));
    
    tmp = (9*gamma2 + 1 - c2)/(2*D2);
    double tmp1 = -9*gamma/D2;
    double a31 = tmp1*(c3*(k*a23 - b21) + k*c4*(1 + 0.25*k2)) + tmp*(3*c3*(2*a23 - k*b21) + c4*(2 + 3*k2));
    double a32 = tmp1*(c3*(k*a24 - b22) + k*0.25*c4) - 3*tmp*(c3*(k*b22 + d21 - 2*a24) - c4);
    
    tmp = (9*gamma2 + 1 + 2*c2)/(8*D2);
    tmp1 = 3*gamma/D2;
    double b31 = tmp1*(3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k2)) + tmp*(12*c3*(k*a23 - b21) + 3*c4*k*(4 + k2));
    double b32 = tmp1*(3*c3*(k*b22 + d21 - 2*a24) - 3*c4) + tmp*(12*c3*(k*a24 - b22) + 3*c4*k);
    
    double a1 = -1.5*c3*(2*a21 - sign*a23 + d21*(2 - 3*sign)) - 0.375*c4*(8 - 4*sign - k2*(2 + sign));
    double a2 = 1.5*c3*(a24 - 2*a22) + 1.125*c4;
    
    tmp = 1/D3;
    double s1 = tmp*(1.5*c3*(2*a21*(k2 - 2) - a23*(k2 + 2) - 2*k*b21) - 0.375*c4*(k2*(3*k2 - 8) + 8) );
    double s2 = tmp*(1.5*c3*(2*a22*(k2 - 2) - sign*a24*(k2 + 2) - 2*sign*k*b22 + d21*(2 - 3*sign)) + 0.375*c4*((8 - 4*sign) - k2*(2 + sign)));
    
    tmp = -k/(16*gamma);
    tmp1 = gamma*(gamma*k - 1);
    double b33 = tmp*(12*c3*(b21 + k*(a23 - 2*a21)) + 3*c4*k*(3*k2 - 4) + 16*s1*tmp1);
    double b34 = 2*tmp*(-12*c3*k*a22 + 3*c4*k + 8*s2*tmp1);
    double b35 = tmp*(12*c3*(b22 + k*a24) + 3*c4*k);
    
    double l1 = a1 + 2*gamma2*s1;
    double l2 = a2 + 2*gamma2*s2;
    
    double delta = gamma2 - c2;
    
    double r1 = yLA*this->sma;		
    double Az = Az0/r1;
    
    double n1 = this->mean_motion; // sqrt(CUB(yLA)); // normalized, (might be a better way to compute this haha)
    double s = n1*time;
    
    double Az2 = Az*Az;
    double Ax2 = -(Az2*l2 + delta)/l1;
    double Ax = sqrt(Ax2);
    
    double omega2 = s1*Ax2 + s2*Az2;
    double omega = 1 + omega2;
    double tau = omega*s;
    
    double tau1 = gamma*tau + phi;
    
    tmp = gamma*omega*n1;
    
    //double period = 6.283185307179586476925286766559/tmp; // correct
    //std::cout << period/86400 << std::endl;
    
    double st1 = sin(tau1);
    double st2 = sin(2*tau1);
    double st3 = sin(3*tau1);
    
    double ct1 = cos(tau1);
    double ct2 = cos(2*tau1);
    double ct3 = cos(3*tau1);
    
    double coef2 = a23*Ax2 + sign*a24*Az2;
    double coef3 = Ax*(a31*Ax2 + sign*a32*Az2);
    double x = a21*Ax2 + a22*Az2 - Ax*ct1 + coef2*ct2 + coef3*ct3  ;
    double u = tmp*( Ax*st1 - coef2*2*st2 - coef3*3*st3 );
    
    double coef1 = k*Ax + Ax*(b33*Ax2 + b34*Az2 + sign*b35*Az2);
    coef2 = b21*Ax2 + sign*b22*Az2;
    coef3 = Ax*(b31*Ax2 + sign*b32*Az2);
    double y = coef1*st1 + coef2*st2 + coef3*st3;
    double v = tmp*( coef1*ct1 + coef2*2*ct2 + coef3*3*ct3 );
    
    double z,w;
    if(n == 1 || n == 3 ) {
        z = Az*ct1 + d21*Ax*Az*(ct2 - 3) + Az*(d32*Ax2 - d31*Az2)*ct3;
        w = -tmp*(Az*st1 + 2*d21*Ax*Az*st2 + 3*Az*(d32*Ax2 - d31*Az2)*st3 );
        if(n == 3) {
            z *= -1;
            w *= -1;
        } 
    } else {
        z = Az*st1 + d21*Ax*Az*st2 + Az*(d32*Ax2 + d31*Az2)*st3;
        w = tmp*(Az*ct1 + 2*d21*Ax*Az*ct2 + 3*Az*(d32*Ax2 + d31*Az2)*ct3 );
        if(n == 2) {
            z *= -1;
            w *= -1;
        }
    }
    
    double ydot = xL1 - (xL1*n1 - yLA*v)/n1;
    double xdot = yLA*u/n1;
    double zdot = yLA*w/n1;
    std::array<double,6> x0 = {x*yLA + xL1,y*yLA,z*yLA,xdot,ydot,zdot};
    return x0;
}

std::array<double,6> CR3BP::convert_state_to_inertial(const std::array<double,6>& x) {
    std::array<double,6> inertial;
    inertial[0] = x[0]*this->sma;
    inertial[1] = x[1]*this->sma;
    inertial[2] = x[2]*this->sma;
    double v = this->sma*this->mean_motion;
    inertial[3] = x[3]*v - inertial[1]*this->mean_motion;
    inertial[4] = x[4]*v + inertial[0]*this->mean_motion;
    inertial[5] = x[5]*v;
    return inertial;
}

std::array<double,6> CR3BP::get_state_rate(const std::array<double,6>& x, const double& t) const {
    double dx1 = x[0] + mu;
    double dx2 = x[0] - mu1;
    double d = x[1]*x[1] + x[2]*x[2];
    double r1_sq = dx1*dx1 + d;
    double r2_sq = dx2*dx2 + d;
    
    double r1 = sqrt(r1_sq);
    double r2 = sqrt(r2_sq);
    
    double g1 = mu1/(r1*r1_sq);
    double g2 = mu/(r2*r2_sq);
    double g = g1 + g2;
    
    double x_dd = x[0] + 2*x[4] - dx1*g1 - dx2*g2;
    double y_dd = x[1]*(1 - g) - 2*x[3];
    double z_dd = -x[2]*g;
    std::array<double,6> dx = {x[3], x[4], x[5], x_dd, y_dd,z_dd};
    return dx;
}

VectorTable::VectorTable(const std::string& fn) {
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
    
    uint_fast32_t nLines = lines.size();
    uint_fast32_t i = 0;
    double x,y,z;
    while(i < nLines) {
        this->jds.push_back(std::stod(lines[i]));
        i++;
        x = std::stod(lines[i].substr(0,23));
        y = std::stod(lines[i].substr(23,46));
        z = std::stod(lines[i].substr(46,69));
        this->pos.emplace_back(std::array<double,3>{x,y,z});
        i++;
        x = std::stod(lines[i].substr(0,23));
        y = std::stod(lines[i].substr(23,46));
        z = std::stod(lines[i].substr(46,69));
        this->vel.emplace_back(std::array<double,3>{x,y,z});
        i++;
    }
}

std::array<double,3> VectorTable::getPos(const double& jd) const {
    if(jd < this->jds[0]){
        return this->pos[0];
    }

    if(jd > this->jds.back()){
        return this->pos.back();
    }
    
    uint_fast32_t idx = Math::lowerbound(this->jds.data(),jd,this->jds.size());
    uint_fast32_t idx1 = idx + 1;
    
    double dt = (jd - this->jds[idx])*86400;
    double dt_hi = (this->jds[idx1] - this->jds[idx])*86400;
    
    double b = dt_hi*dt_hi;
    double a = b*dt_hi;
    double c = 3*b;
    double d = 2*dt_hi;
    double det = 1/(a*d - b*c);
    
    const std::array<double,3>& pos_lo = this->pos[idx];
    const std::array<double,3>& pos_hi = this->pos[idx1];
    const std::array<double,3>& vel_lo = this->vel[idx];
    const std::array<double,3>& vel_hi = this->vel[idx1];
    
    std::array<double,3> pos;
    for(int i = 0; i < 3; i++){
        double x = pos_hi[i] - pos_lo[i] - vel_lo[i]*dt_hi;
        double v = vel_hi[i] - vel_lo[i];
        double coef1 = d*x - b*v;
        double coef2 = a*v - c*x;
        pos[i] = pos_lo[i] + dt*(vel_lo[i] + dt*det*(coef2 + dt*coef1));
    }
    return pos;
}

std::array<double,3> VectorTable::getVel(const double& jd) const {
    if(jd < this->jds[0]){
        return this->vel[0];
    }

    if(jd > this->jds.back()){
        return this->vel.back();
    }
    
    uint_fast32_t idx = Math::lowerbound(this->jds.data(),jd,this->jds.size());
    
    double f_hi = (jd - this->jds[idx])/(this->jds[idx + 1] - this->jds[idx]);
    double f_lo = 1.0 - f_hi;
    
    const std::array<double,3>& vel_lo = this->vel[idx];
    const std::array<double,3>& vel_hi = this->vel[idx + 1];
    std::array<double,3> vel;
    for(int i = 0; i < 3; i++){
        vel[i] = vel_lo[i]*f_lo + vel_hi[i]*f_hi;
    }
    return vel;
}

EarthMoonSun::EarthMoonSun(const double& JD0) : jd0(JD0) {
    this->earth = new VectorTable("resources/EARTH_EMB_VECTOR.dat");
    this->moon = new VectorTable("resources/MOON_EMB_VECTOR.dat");
    this->emb = new VectorTable("resources/EMB_SSB_VECTOR.dat");
    this->sun = new VectorTable("resources/SUN_SSB_VECTOR.dat");
    
    this->jdf = this->earth->jds.back();
    
    double jd = jd0 - 2;
    //int nTable = (jdf - jd0)/djd;
    while(jd < jdf){
        this->jds.push_back(jd);
        this->r_earth.push_back(this->earth->getPos(jd));
        this->r_moon.push_back(this->moon->getPos(jd));
        std::array<double,3> emb_ssb = this->emb->getPos(jd);
        this->r_emb.push_back(emb_ssb);
        std::array<double,3> sun_ssb = this->sun->getPos(jd);
        std::array<double,3> sun_emb = {sun_ssb[0] - emb_ssb[0],sun_ssb[1] - emb_ssb[1],sun_ssb[2] - emb_ssb[2] }; // GCRS and ICRF close enough
        this->r_sun.push_back(sun_emb);
        
        std::array<double,3> emb_ssb_v = this->emb->getVel(jd);
        this->v_emb.push_back(emb_ssb_v);
        
        std::array<double,3> h = Math::cross(emb_ssb,emb_ssb_v);
        double r2 = Math::dot(emb_ssb,emb_ssb);
        h[0] /= r2;
        h[1] /= r2;
        h[2] /= r2;
        
        this->omega.push_back(h);
        jd += djd;
    }
    
    int n = this->jds.size() - 1;
    for(int i = 0; i < n;i++) {
        const std::array<double,3>& p1 = this->r_earth[i];
        const std::array<double,3>& p2 = this->r_earth[i+1];
        std::array<double,3> dp;
        for(int j = 0; j < 3;j++) {
            dp[j] = (p2[j] - p1[j])/djd;
        }
        this->d_earth.push_back(dp);
    }
    
    for(int i = 0; i < n;i++) {
        const std::array<double,3>& p1 = this->r_moon[i];
        const std::array<double,3>& p2 = this->r_moon[i+1];
        std::array<double,3> dp;
        for(int j = 0; j < 3;j++) {
            dp[j] = (p2[j] - p1[j])/djd;
        }
        this->d_moon.push_back(dp);
    }
    
    for(int i = 0; i < n;i++) {
        const std::array<double,3>& p1 = this->r_sun[i];
        const std::array<double,3>& p2 = this->r_sun[i+1];
        std::array<double,3> dp;
        for(int j = 0; j < 3;j++) {
            dp[j] = (p2[j] - p1[j])/djd;
        }
        this->d_sun.push_back(dp);
    }
    
    for(int i = 0; i < n;i++) {
        const std::array<double,3>& p1 = this->r_emb[i];
        const std::array<double,3>& p2 = this->r_emb[i+1];
        std::array<double,3> dp;
        for(int j = 0; j < 3;j++) {
            dp[j] = (p2[j] - p1[j])/djd;
        }
        this->d_emb.push_back(dp);
    }
    
    for(int i = 0; i < n;i++) {
        const std::array<double,3>& p1 = this->v_emb[i];
        const std::array<double,3>& p2 = this->v_emb[i+1];
        std::array<double,3> dp;
        for(int j = 0; j < 3;j++) {
            dp[j] = (p2[j] - p1[j])/djd;
        }
        this->dv_emb.push_back(dp);
    }
}

std::array< std::array<double,3>, 4> EarthMoonSun::getEarthMoonBarycenterCS(const double& jd) const {
    std::array<double,3> r_moon =  this->moon->getPos(jd);
    std::array<double,3> r_earth =  this->earth->getPos(jd);
    std::array<double,3> x = {r_moon[0] - r_earth[0],r_moon[1] - r_earth[1],r_moon[2] - r_earth[2]};
    //Normalized: Earth Distance from barycetner
    double xL1 = OrbitalElements::MOON_MU/(OrbitalElements::EARTH_MU + OrbitalElements::MOON_MU);
    std::array<double,3> Barycenter = {r_earth[0] + xL1*x[0],r_earth[1] + xL1*x[1],r_earth[2] + xL1*x[2]};
    
    std::array<double,3> v_earth = this->earth->getVel(jd);
    std::array<double,3> h_earth = Math::cross(r_earth,v_earth);
    std::array<double,3> v_moon = this->moon->getVel(jd);
    std::array<double,3> h_moon = Math::cross(r_moon,v_moon);
    
    std::array<double,3> z = {h_earth[0] + h_moon[0],h_earth[1] + h_moon[1],h_earth[2] + h_moon[2]};
    double z_mag = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    double x_mag = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    for(int i = 0;i < 3; i++) {
        x[i] /= x_mag;
        z[i] /= z_mag;
    }
    std::array<double,3> y = Math::cross(z,x);
    std::array< std::array<double,3>, 4> CS = {Barycenter,x,y,z};
    return CS;
}

std::array< std::array<double,3>, 4> EarthMoonSun::getEarthMoonL1CS(const double& jd) const {
    std::array<double,3> r_moon =  this->moon->getPos(jd);
    std::array<double,3> r_earth =  this->earth->getPos(jd);
    std::array<double,3> x = {r_moon[0] - r_earth[0],r_moon[1] - r_earth[1],r_moon[2] - r_earth[2]};
    //Normalized: L1 from  Earth Distance
    
    double xL1 = OrbitalElements::MOON_MU/(OrbitalElements::EARTH_MU + OrbitalElements::MOON_MU) + CR3BP::EARTH_L1;

    std::array<double,3> L1 = {r_earth[0] + xL1*x[0],r_earth[1] + xL1*x[1],r_earth[2] + xL1*x[2]};
    
    std::array<double,3> v_earth = this->earth->getVel(jd);
    std::array<double,3> h_earth = Math::cross(r_earth,v_earth);
    std::array<double,3> v_moon = this->moon->getVel(jd);
    std::array<double,3> h_moon = Math::cross(r_moon,v_moon);
    
    std::array<double,3> z = {h_earth[0] + h_moon[0],h_earth[1] + h_moon[1],h_earth[2] + h_moon[2]};
    double z_mag = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    double x_mag = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    for(int i = 0;i < 3; i++) {
        x[i] /= x_mag;
        z[i] /= z_mag;
    }
    std::array<double,3> y = Math::cross(z,x);
    std::array< std::array<double,3>, 4> CS = {L1,x,y,z};
    return CS;
}

std::array<double,6> EarthMoonSun::cr3bp_to_embj2000(const std::array<double,6>& cr3bp, const double& jd) const {
    std::array< std::array<double,3>, 4> frame = this->getEarthMoonBarycenterCS(jd);
    std::array< std::array<double,3>, 3> CS = {frame[1],frame[2],frame[3]};
    std::array< std::array<double,3>, 3> CST = Math::transpose(CS);
    std::array<double,3> origin = frame[0];
    std::array<double,3> pos = {origin[0] + cr3bp[0],origin[1] + cr3bp[1],origin[2] + cr3bp[2]};
    std::array<double,3> vel = {cr3bp[3],cr3bp[4],cr3bp[5]};
    
    pos = Math::mult(CST,pos);
    vel = Math::mult(CST,vel);

    const std::array<double,6> x = {pos[0],pos[1],pos[2],vel[0],vel[1],vel[2]};
    return x;
}

void EarthMoonSun::getA(const std::array<double,6>& x, const double& jd, double** A) const {
    Math::clear(A,6);
    A[0][3] = 1;
    A[1][4] = 1;
    A[2][5] = 1;
    
    std::size_t idx = (std::size_t)((jd - this->jd0)/this->djd);
    if(idx < 0){
        idx = 0;
    }
    if(idx > this->jds.size() - 2) {
        std::cout << jd << std::endl;
        throw "bad jd";
        // could compute from table if bad
    }
    double djd = jd - this->jds[idx];
    
    const std::array<double,3>& earth = this->r_earth[idx];
    const std::array<double,3>& moon = this->r_moon[idx];
    const std::array<double,3>& sun = this->r_sun[idx];
    const std::array<double,3>& emb = this->r_emb[idx];
    //const std::array<double,3>& vemb = this->v_emb[idx];
    
    const std::array<double,3>& dE = this->d_earth[idx];
    const std::array<double,3>& dM = this->d_moon[idx];
    const std::array<double,3>& dS = this->d_sun[idx];
    const std::array<double,3>& dEMB = this->d_emb[idx];
    //const std::array<double,3>& dvEMB = this->dv_emb[idx];
    
    std::array<double,3> rE;
    std::array<double,3> rM;
    std::array<double,3> rS;
    std::array<double,3> rEMB;
    //std::array<double,3> vEMB;
    
    for(int i = 0; i < 3; i++){
        rE[i] = x[i] - (earth[i] + dE[i]*djd);
        rM[i] = x[i] - (moon[i] + dM[i]*djd);
        rS[i] = x[i] - (sun[i] + dS[i]*djd);
        rEMB[i] = emb[i] + dEMB[i]*djd;
        // vEMB[i] = 2*(vemb[i] + dvEMB[i]*djd); // pre mult
    }

    double r_E_sq = Math::dot(rE,rE);
    double r_E = sqrt(r_E_sq);
    double g_E = OrbitalElements::EARTH_MU/(r_E*r_E_sq);
    double mu_E = 3*g_E/r_E_sq;
    

    double r_M_sq = Math::dot(rM,rM);
    double r_M = sqrt(r_M_sq);
    double g_M = OrbitalElements::MOON_MU/(r_M*r_M_sq);
    double mu_M = 3*g_M/r_M_sq;

    double r_S_sq = Math::dot(rS,rS);
    double r_S = sqrt(r_S_sq);
    double g_S = OrbitalElements::SUN_MU/(r_S*r_S_sq);
    double mu_S = 3*g_S/r_S_sq;

    double g = g_E + g_M + g_S;
    
    const std::array<double,3>& w = this->omega[idx]; 
    // cross : w[1]*v[2] - w[2]*v[1], w[2]*v[0] - w[0]*v[2], w[0]*v[1] - w[1]*v[0];
    // dr/dx = {1,0,0} omega x drdx = {0,w[2],-w[1]} = v omega x omega x drdx = 0,w[0]*w[1],w[0]*w[2] ie. {0,wxy,wxz}
    // dr/dy = {0,1,0} omega x drdy = {-w[2],0,w[0]} = v omega x omega x drdy = w[1]*w[0],0,w[1]*w[2] ie. {wxy,0,wyz};
    // dr/dz = {0,0,1} omega x drdz = {w[1],-w[0],0} = v omega x omega x drdz = w[2]*w[0],w[2]*w[1],0 ie. {wxz,wyz,0};
    double wyz = w[1]*w[2];
    double wxz = w[0]*w[2];
    double wxy = w[0]*w[1];
    
    // dv dx is zero, so coriolos component is zero
    
    A[3][0] = mu_E*rE[0]*rE[0] + mu_M*rM[0]*rM[0] + mu_S*rS[0]*rS[0] - g ;
    A[3][1] = mu_E*rE[0]*rE[1] + mu_M*rM[0]*rM[1] + mu_S*rS[0]*rS[1] + wxy;
    A[3][2] = mu_E*rE[0]*rE[2] + mu_M*rM[0]*rM[2] + mu_S*rS[0]*rS[2] + wxz;
    
    A[4][0] = A[3][1];
    A[4][1] = mu_E*rE[1]*rE[1] + mu_M*rM[1]*rM[1] + mu_S*rS[1]*rS[1] - g ;
    A[4][2] = mu_E*rE[1]*rE[2] + mu_M*rM[1]*rM[2] + mu_S*rS[1]*rS[2] + wyz;
    
    A[5][0] = A[3][2];
    A[5][1] = A[4][2];
    A[5][2] = mu_E*rE[2]*rE[2] + mu_M*rM[2]*rM[2] + mu_S*rS[2]*rS[2] - g;
}

void EarthMoonSun::getA(const std::array<double,6>& x, const double& jd, Matrix<6,6>& A) const {
    A.set_zero();
    A[0][3] = 1;
    A[1][4] = 1;
    A[2][5] = 1;
    
    std::size_t idx = (std::size_t)((jd - this->jd0)/this->djd);
    if(idx < 0){
        idx = 0;
    }
    if(idx > this->jds.size() - 2) {
        std::cout << jd << std::endl;
        throw "bad jd";
        // could compute from table if bad
    }
    double djd = jd - this->jds[idx];
    
    const std::array<double,3>& earth = this->r_earth[idx];
    const std::array<double,3>& moon = this->r_moon[idx];
    const std::array<double,3>& sun = this->r_sun[idx];
    const std::array<double,3>& emb = this->r_emb[idx];
    //const std::array<double,3>& vemb = this->v_emb[idx];
    
    const std::array<double,3>& dE = this->d_earth[idx];
    const std::array<double,3>& dM = this->d_moon[idx];
    const std::array<double,3>& dS = this->d_sun[idx];
    const std::array<double,3>& dEMB = this->d_emb[idx];
    //const std::array<double,3>& dvEMB = this->dv_emb[idx];
    
    std::array<double,3> rE;
    std::array<double,3> rM;
    std::array<double,3> rS;
    std::array<double,3> rEMB;
    //std::array<double,3> vEMB;
    
    for(int i = 0; i < 3; i++){
        rE[i] = x[i] - (earth[i] + dE[i]*djd);
        rM[i] = x[i] - (moon[i] + dM[i]*djd);
        rS[i] = x[i] - (sun[i] + dS[i]*djd);
        rEMB[i] = emb[i] + dEMB[i]*djd;
        // vEMB[i] = 2*(vemb[i] + dvEMB[i]*djd); // pre mult
    }

    double r_E_sq = Math::dot(rE,rE);
    double r_E = sqrt(r_E_sq);
    double g_E = OrbitalElements::EARTH_MU/(r_E*r_E_sq);
    double mu_E = 3*g_E/r_E_sq;
    

    double r_M_sq = Math::dot(rM,rM);
    double r_M = sqrt(r_M_sq);
    double g_M = OrbitalElements::MOON_MU/(r_M*r_M_sq);
    double mu_M = 3*g_M/r_M_sq;

    double r_S_sq = Math::dot(rS,rS);
    double r_S = sqrt(r_S_sq);
    double g_S = OrbitalElements::SUN_MU/(r_S*r_S_sq);
    double mu_S = 3*g_S/r_S_sq;

    double g = g_E + g_M + g_S;
    
    const std::array<double,3>& w = this->omega[idx]; 
    // cross : w[1]*v[2] - w[2]*v[1], w[2]*v[0] - w[0]*v[2], w[0]*v[1] - w[1]*v[0];
    // dr/dx = {1,0,0} omega x drdx = {0,w[2],-w[1]} = v omega x omega x drdx = 0,w[0]*w[1],w[0]*w[2] ie. {0,wxy,wxz}
    // dr/dy = {0,1,0} omega x drdy = {-w[2],0,w[0]} = v omega x omega x drdy = w[1]*w[0],0,w[1]*w[2] ie. {wxy,0,wyz};
    // dr/dz = {0,0,1} omega x drdz = {w[1],-w[0],0} = v omega x omega x drdz = w[2]*w[0],w[2]*w[1],0 ie. {wxz,wyz,0};
    double wyz = w[1]*w[2];
    double wxz = w[0]*w[2];
    double wxy = w[0]*w[1];
    
    // dv dx is zero, so coriolos component is zero
    
    A[3][0] = mu_E*rE[0]*rE[0] + mu_M*rM[0]*rM[0] + mu_S*rS[0]*rS[0] - g ;
    A[3][1] = mu_E*rE[0]*rE[1] + mu_M*rM[0]*rM[1] + mu_S*rS[0]*rS[1] + wxy;
    A[3][2] = mu_E*rE[0]*rE[2] + mu_M*rM[0]*rM[2] + mu_S*rS[0]*rS[2] + wxz;
    
    A[4][0] = A[3][1];
    A[4][1] = mu_E*rE[1]*rE[1] + mu_M*rM[1]*rM[1] + mu_S*rS[1]*rS[1] - g ;
    A[4][2] = mu_E*rE[1]*rE[2] + mu_M*rM[1]*rM[2] + mu_S*rS[1]*rS[2] + wyz;
    
    A[5][0] = A[3][2];
    A[5][1] = A[4][2];
    A[5][2] = mu_E*rE[2]*rE[2] + mu_M*rM[2]*rM[2] + mu_S*rS[2]*rS[2] - g;
}

Matrix<3,3> EarthMoonSun::getG(const std::array<double,6>& x, const double& jd) const {

    std::size_t idx = (std::size_t)((jd - this->jd0)/this->djd);
    if(idx < 0){
        idx = 0;
    }
    if(idx > this->jds.size() - 2) {
        idx = this->jds.size() - 2;
    }
    double djd = jd - this->jds[idx];
    
    const std::array<double,3>& earth = this->r_earth[idx];
    const std::array<double,3>& moon = this->r_moon[idx];
    const std::array<double,3>& sun = this->r_sun[idx];
    const std::array<double,3>& emb = this->r_emb[idx];
    
    const std::array<double,3>& dE = this->d_earth[idx];
    const std::array<double,3>& dM = this->d_moon[idx];
    const std::array<double,3>& dS = this->d_sun[idx];
    const std::array<double,3>& dEMB = this->d_emb[idx];
    
    std::array<double,3> rE;
    std::array<double,3> rM;
    std::array<double,3> rS;
    std::array<double,3> rEMB;
    
    for(int i = 0; i < 3; i++){
        rE[i] = x[i] - (earth[i] + dE[i]*djd);
        rM[i] = x[i] - (moon[i] + dM[i]*djd);
        rS[i] = x[i] - (sun[i] + dS[i]*djd);
        rEMB[i] = emb[i] + dEMB[i]*djd;
    }

    double r_E_sq = Math::dot(rE,rE);
    double r_E = sqrt(r_E_sq);
    double g_E = OrbitalElements::EARTH_MU/(r_E*r_E_sq);
    double mu_E = 3*g_E/r_E_sq;

    double r_M_sq = Math::dot(rM,rM);
    double r_M = sqrt(r_M_sq);
    double g_M = OrbitalElements::MOON_MU/(r_M*r_M_sq);
    double mu_M = 3*g_M/r_M_sq;

    double r_S_sq = Math::dot(rS,rS);
    double r_S = sqrt(r_S_sq);
    double g_S = OrbitalElements::SUN_MU/(r_S*r_S_sq);
    double mu_S = 3*g_S/r_S_sq;

    double g = g_E + g_M + g_S;
    
    const std::array<double,3>& w = this->omega[idx]; 
    double wyz = w[1]*w[2];
    double wxz = w[0]*w[2];
    double wxy = w[0]*w[1];
    Matrix<3,3> G;
    G.data[0] = mu_E*rE[0]*rE[0] + mu_M*rM[0]*rM[0] + mu_S*rS[0]*rS[0] - g ;
    G.data[1] = mu_E*rE[0]*rE[1] + mu_M*rM[0]*rM[1] + mu_S*rS[0]*rS[1] + wxy;
    G.data[2] = mu_E*rE[0]*rE[2] + mu_M*rM[0]*rM[2] + mu_S*rS[0]*rS[2] + wxz;
    
    G.data[3] = G.data[1];
    G.data[4] = mu_E*rE[1]*rE[1] + mu_M*rM[1]*rM[1] + mu_S*rS[1]*rS[1] - g ;
    G.data[5] = mu_E*rE[1]*rE[2] + mu_M*rM[1]*rM[2] + mu_S*rS[1]*rS[2] + wyz;
    
    G.data[6] = G.data[2];
    G.data[7] = G.data[5];
    G.data[8] = mu_E*rE[2]*rE[2] + mu_M*rM[2]*rM[2] + mu_S*rS[2]*rS[2] - g;
    return G;
}

std::array<double,6> EarthMoonSun::get_state_rate(const std::array<double,6>& x, const double& t) const {
    double jd = this->jd0 + t/86400.0;
    
    uint_fast16_t idx = (uint_fast16_t)((jd - this->jd0)/this->djd);
    if(idx < 0){
        idx = 0;
    }
    if(idx > this->jds.size() - 2) {
        std::cout << jd << std::endl;
        throw "bad jd";
    }
    double djd = jd - this->jds[idx];
    
    const std::array<double,3>& earth = this->r_earth[idx];
    const std::array<double,3>& moon = this->r_moon[idx];
    const std::array<double,3>& sun = this->r_sun[idx];
    const std::array<double,3>& emb = this->r_emb[idx];
    //const std::array<double,3>& vemb = this->v_emb[idx];
    
    const std::array<double,3>& dE = this->d_earth[idx];
    const std::array<double,3>& dM = this->d_moon[idx];
    const std::array<double,3>& dS = this->d_sun[idx];
    const std::array<double,3>& dEMB = this->d_emb[idx];
    //const std::array<double,3>& dvEMB = this->dv_emb[idx];
    
    std::array<double,3> rE;
    std::array<double,3> rM;
    std::array<double,3> rS;
    std::array<double,3> rEMB;
    //std::array<double,3> vEMB;
    
    for(uint_fast8_t i = 0; i < 3; i++){
        rE[i] = x[i] - (earth[i] + dE[i]*djd);
        rM[i] = x[i] - (moon[i] + dM[i]*djd);
        rS[i] = x[i] - (sun[i] + dS[i]*djd);
        rEMB[i] = emb[i] + dEMB[i]*djd;
        // vEMB[i] = 2*(vemb[i] + dvEMB[i]*djd); // pre mult
    }
    
    double r_E_sq = Math::dot(rE,rE);
    double g_E = OrbitalElements::EARTH_MU/(sqrt(r_E_sq)*r_E_sq);
    
    double r_M_sq = Math::dot(rM,rM);
    double g_M = OrbitalElements::MOON_MU/(sqrt(r_M_sq)*r_M_sq);
    
    double r_S_sq = Math::dot(rS,rS);
    double g_S = OrbitalElements::SUN_MU/(sqrt(r_S_sq)*r_S_sq);
    
    const std::array<double,3>& w = this->omega[idx]; 
    
    const std::array<double,3>& a_cen = Math::cross(w,Math::cross(w,rEMB));
    
    std::array<double,3> v = {x[3], x[4], x[5]};
    const std::array<double,3>& a_cor = Math::cross(w,v); // by definition coriolos effect from frame should be zero since origin of frame (barycenter) determines angular velocity (mean motion)
    
    double x_dd = -g_E*rE[0] - g_M*rM[0] - g_S*rS[0] - a_cen[0] - 2*a_cor[0];
    double y_dd = -g_E*rE[1] - g_M*rM[1] - g_S*rS[1] - a_cen[1] - 2*a_cor[1];
    double z_dd = -g_E*rE[2] - g_M*rM[2] - g_S*rS[2] - a_cen[2] - 2*a_cor[2];
    std::array<double,6> dx = {x[3], x[4], x[5], x_dd, y_dd,z_dd};
    return dx;
}
