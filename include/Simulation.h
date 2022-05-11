#pragma once

#include <array>
#include <vector>
#include "Matrix.h"
#include "Dynamics.h"
#include "Ode.h"
#include <thread>

struct Section {
	std::array<double,6> initial_state;
	double t_start, t_final;
	std::vector<double> times;
	std::vector< std::array<double,6> > states;
	Matrix<6,6> STM;
	ODE_RK4<6> ode;
	EarthMoonSun* dynamics;
	static const double SECTION_DAYS;
	
	Section(EarthMoonSun* f) {
		this->dynamics = f;
		this->ode = ODE_RK4<6>();
		this->ode.set_dynamics(f);
		this->ode.recording.set_record_interval(100);
		this->ode.set_timestep(10);
	}
	
	void compute_states(){
		this->ode.recording.clear();
		this->ode.set_time(this->t_start);
		this->ode.run(this->initial_state,this->t_final);
		this->times = this->ode.recording.get_times();
		this->states = this->ode.recording.get_states();
		this->times.push_back(ode.get_time());
		this->states.push_back(ode.get_state());
	}

    void target(const std::array<double,6>& xp);

	void compute_STM();

    static void minimizeDX(std::vector<Section>& sections);

    static void minimizeDV(std::vector<Section>& sections);

    static double calcDV(std::vector<Section>& sections);
};

const double Section::SECTION_DAYS = 3;

std::array<double,6> convert(CR3BP* cr3bp, EarthMoonSun* dynamics, const std::array<double,6>& state_guess, const double& jd){
	std::array<double,6> x = cr3bp->convert_state_to_inertial(state_guess);
	
	const double xL1 = cr3bp->getL1() - cr3bp->mu;
	
	std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonBarycenterCS(jd);
	
	std::array<double,3> moon = dynamics->moon->getPos(jd - 0.0005);
	std::array<double,3> earth = dynamics->earth->getPos(jd - 0.0005);
	std::array<double,3> L1_0 = {moon[0] - earth[0],moon[1] - earth[1],moon[2] - earth[2]};
	
	moon = dynamics->moon->getPos(jd + 0.0005);
	earth = dynamics->earth->getPos(jd + 0.0005);
	std::array<double,3> L1_1 = {moon[0] - earth[0],moon[1] - earth[1],moon[2] - earth[2]};
	
	double dt = 0.001*86400;
	double dL1dt = (sqrt(Math::dot(L1_1,L1_1)) - sqrt(Math::dot(L1_0,L1_0)))*xL1/dt;
	
	//std::cout << jd << ": " << dL1 << std::endl;
	
	std::array< std::array<double,3>, 3> CS = {frame[1],frame[2],frame[3]};
	std::array< std::array<double,3>, 3> CST = Math::transpose(CS);
	
	std::array<double,3> pos = {x[0],x[1],x[2]};
	std::array<double,3> vel = {x[3] + dL1dt,x[4],x[5]};
	
	pos = Math::mult(CST,pos);
	vel = Math::mult(CST,vel);

	std::array<double,3>& origin = frame[0]; // Should be zero!
	x[0] = pos[0] + origin[0];
	x[1] = pos[1] + origin[1];
	x[2] = pos[2] + origin[2];
	x[3] = vel[0];
	x[4] = vel[1];
	x[5] = vel[2];
	
	return x;
}

double minimizeDV_Gradient(std::vector<Section>& sections){
	const double deltaX = 1;
	const double deltaX_half = deltaX*0.5;
	const double deltaT = 120;
	const double deltaT_half = deltaT*0.5;
	
	const uint_fast16_t nSections = sections.size();
	std::cout << "Minimizing..." << std::endl;
	
	const uint_fast16_t n1 = nSections-1;
	
	double d;

	Matrix<3,3> dvdx;
	Vector<3> dvdt;
	Matrix<7,7> B;
	Vector<7> x;
	std::array<double,6> xf0, xf1;
	std::array<double,6> xf_last = {0};
	double err = 0;

	int iter = 0;
	while(iter++ < 5) {
		for (uint_fast16_t section = 0; section < n1; section++) {
			Section& section_i = sections[section];
			
			const std::array<double,6>& x_front = sections[section+1].initial_state;

			if(section > 0){
				xf_last = sections[section-1].states.back();
			} else {
				xf_last = section_i.states[0];
			}
				
			section_i.compute_states();
			
			xf1  = section_i.states.back();

			x.set_zero();
			x.data[4] = -(x_front[3] - xf1[3]);
			x.data[5] = -(x_front[4] - xf1[4]);
			x.data[6] = -(x_front[5] - xf1[5]);
			
			d = x.data[4]*x.data[4] + x.data[5]*x.data[5] + x.data[6]*x.data[6];

			if(d < 1e-6) {
				continue;
			}
			
			for(uint_fast8_t i = 0; i < 3; i++){
				section_i.initial_state[i] -= deltaX_half;
				section_i.compute_states();
				
				xf0  = section_i.states.back();
				
				section_i.initial_state[i] += deltaX;
				section_i.compute_states();
				
				xf1  = section_i.states.back();
				
				section_i.initial_state[i] -= deltaX_half;
				
				dvdx[0][i] = (xf1[3] - xf0[3])/deltaX;
				dvdx[1][i] = (xf1[4] - xf0[4])/deltaX;
				dvdx[2][i] = (xf1[5] - xf0[5])/deltaX;
			}
			
			section_i.t_start -= deltaT_half;
			section_i.compute_states();
			
			xf0 = section_i.states.back();

			section_i.t_start += deltaT;
			section_i.compute_states();
			
			xf1  = section_i.states.back();
			
			section_i.t_start -= deltaT_half;
			
			dvdt[0] = (xf1[3] - xf0[3])/deltaT;
			dvdt[1] = (xf1[4] - xf0[4])/deltaT;
			dvdt[2] = (xf1[5] - xf0[5])/deltaT;

			// we should change initial state as little as possible so
			// minimize dl = (vx_last*dt - dx)^2 + (vy_last*dt - dy)^2 + (vz_last*dt - dz)^2
			//  dvx + dx*dvxdx + dy*dvxdy + dz*dvxdz + dt*dvxdt = 0
			//  dvy + dx*dvydx + dy*dvydy + dz*dvydz + dt*dvydt = 0
			//  dvz + dx*dvzdx + dy*dvzdy + dz*dvzdz + dt*dvzdt = 0
			// x vector = [gamma1 gamma2 gamma3 dx dy dz dt ] to prevent zeros in diagonal
			B.set_zero();
			double* row;
			for(uint16_t i = 0; i < 3; i++) {
				row = B[i];
				row[0] = dvdx[0][i];
				row[1] = dvdx[1][i];
				row[2] = dvdx[2][i];
				row[3 + i] = 2;
				row[6] = -2*xf_last[3 + i];
			}
			
			row = B[3];
			row[6] = 0;
			for(uint16_t i = 0; i < 3; i++) {
				row[i] = dvdt[i];
				uint16_t idx = 3+i;
				row[idx] = -2*xf_last[idx];
				row[6] -= row[idx]*xf_last[idx];
			}

			for(uint16_t i = 0; i < 3; i++) {
				row = B[4+i];
				row[3] = dvdx[i][0];
				row[4] = dvdx[i][1];
				row[5] = dvdx[i][2];
				row[6] = dvdt[i];
			}
			
			// 
			try {
				Matrix<7,7>::solve(B,x);
			} catch (...) {
				std::cout << "Singular Matrix, simplifying problem" << std::endl;

				xf1  = section_i.states.back();
				Vector<3> dv;
				dv.data[0] = -(x_front[3] - xf1[3]);
				dv.data[1] = -(x_front[4] - xf1[4]);
				dv.data[2] = -(x_front[5] - xf1[5]);

				Vector<3> dx = MatrixX::invert(dvdx)*dv;

				x[3] = dx[0];
				x[4] = dx[1];
				x[5] = dx[2];
				x[6] = 0;
			}

			for(uint_fast8_t i = 0; i < 3; i++){
				section_i.initial_state[i] += x[3+i];
			}
			section_i.t_start += x[6];

			for(uint_fast8_t i = 0; i < 3; i++){
				std::cout << section_i.initial_state[i] << " ";
			}
			std::cout << section_i.t_start;
			std::cout << std::endl;
			err += d;
		}
	}
	
	for (uint_fast16_t section = 0; section < n1; section++) {
		sections[section].t_final = sections[section+1].t_start;
	}
	
	return err;
}
