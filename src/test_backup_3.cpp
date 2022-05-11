#include "Dynamics.h"
#include "ODE.h"
#include "Math.h"
#include "Util.h"
#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
//#include "../Eigen/Dense"

//typedef Eigen::Matrix<double,6,6> Matrix6d;
//typedef Eigen::Matrix<double,6,1> Vector6d;

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
		this->ode.recording.set_record_interval(150);
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
	
	void compute_STM(){		
		Matrix<6,6> A;
		Matrix<6,6> dSTM;
		this->STM.set_identity();
		A[0][3] = 1;
		A[1][4] = 1;
		A[2][5] = 1;
		
		const std::size_t n = this->times.size() - 1;
		std::array<double,6> a;
		std::array<double,6> x;
		std::array<double,6> ax;
		std::array<double,6> ay;
		std::array<double,6> az;
        for(std::size_t i = 0; i < n;i++) {
			double dt = this->times[i+1] - this->times[i];
			x = this->states[i];
			a = this->dynamics->get_state_rate(x, this->times[i]);
			x[0] += 1e-3;
			ax = this->dynamics->get_state_rate(x, this->times[i]);
			x[0] -= 1e-3;
			x[1] += 1e-3;
			ay = this->dynamics->get_state_rate(x, this->times[i]);
			x[1] -= 1e-3;
			x[2] += 1e-3;
			az = this->dynamics->get_state_rate(x, this->times[i]);
			
			for(std::size_t j = 3; j < 6; j++){
				ax[j] = (ax[j] - a[j])/1e-3;
				ay[j] = (ay[j] - a[j])/1e-3;
				az[j] = (az[j] - a[j])/1e-3;
			}
			
			for(std::size_t j = 0; j < 3; j++){
				std::size_t j3 = j + 3;
				A[3][j] = ax[j3];
				A[4][j] = ay[j3];
				A[5][j] = az[j3];
			}
			
			dSTM = A*STM;
			STM += (A*STM)*dt;
		}
	}
	
	void compute_STM2(){
		Matrix<6,6> A;
		this->STM.set_identity();
		const int n = this->times.size()-1;
        for(int i = 0; i < n;i++) {
			double dt = this->times[i+1] - this->times[i];
			double jd = this->dynamics->getJD0() + this->times[i]/86400;
            this->dynamics->getA(this->states[i],jd,A);
			STM += A*STM*dt;
		}
	}
	
	void compute_STM3(){		
		// Markley Method
		Matrix<6,6> dSTM;
		Matrix<3,3> RR;
		Matrix<3,3> RV;
		Matrix<3,3> VR;
		Matrix<3,3> VV;

		double I[9] = {1,0,0,0,1,0,0,0,1};
		this->STM.set_identity();
		const int n = this->times.size() - 1;
		std::vector< Matrix<3,3> > Gs;
		for (int i = 0; i <= n;i++){
			double jd = this->dynamics->getJD0() + this->times[i]/86400;
			Gs.push_back(std::move(this->dynamics->getG(this->states[i],jd)));
		}	

		for (int i = 0; i < n;i++) {
			double dt = this->times[i+1] - this->times[i];
			double dt_half = dt*0.5;
			double dt6 = dt*dt*0.1666666666666666666666;
			Matrix<3,3>& G0 = Gs[i];
			Matrix<3,3>& G = Gs[i + 1];
			for (int j = 0; j < 9; j++) {
				double tmp = G0.data[j] + G.data[j];
				VR.data[j] = tmp*dt_half;
				RR.data[j] = I[j] + (tmp + G0.data[j])*dt6;
				RV.data[j] = I[j]*dt + VR.data[j]*dt6;
				VV.data[j] = I[j] + (tmp + G.data[j])*dt6; 
			}
			
			for (int j = 3; j < 3; j++) {
				for (int k = 3; k < 3; k++) {
					dSTM.set(i,j,RR(i,j));
					dSTM.set(i,j+3,RV(i,j));
					dSTM.set(i+3,j,VR(i,j));
					dSTM.set(i+3,j+3,VV(i,j));
				}
			}
			
			this->STM = dSTM*this->STM;
		}
	}
};

const double Section::SECTION_DAYS = 2;

void printOut(const Recording<6>&  record, std::string filename){
	FILE * pFile;

	pFile = fopen (filename.c_str(),"w");
	const int n = record.number_entries();
	std::string formatStr = "%12.6f %16.14f %16.14f %16.14f\n";
	const char* format = formatStr.c_str();
	for (int i = 0 ; i < n; i++) {
		const std::array<double,6>& state = record.get(i);
		fprintf (pFile, format,record.time_at(i), state[0], state[1], state[2]);
	}
	fclose (pFile);
}

void printOut(const std::vector<double>& t,const std::vector<std::array<double,3> >& x, std::string filename){
	FILE * pFile;

	pFile = fopen (filename.c_str(),"w");
	const int n = t.size();
	std::string formatStr = "%12.6f %16.14f %16.14f %16.14f\n";
	const char* format = formatStr.c_str();
	for (int i = 0 ; i < n; i++) {
		fprintf (pFile, format,t[i], x[i][0], x[i][1], x[i][2]);
	}
	fclose (pFile);
}

void printOut(EarthMoonSun* dynamics, const std::vector<Section>& sections, std::string filename){
	FILE * pFile;

	pFile = fopen (filename.c_str(),"w");

	std::string formatStr = "%12.6f %16.14f %16.14f %16.14f %16.14f %16.14f %16.14f\n";
	const char* format = formatStr.c_str();
	
	for(const Section& section : sections) {
		const int n = section.times.size();
		for (int i = 0 ; i < n; i++) {
			const std::array<double,6>& x = section.states[i];
			double jd = dynamics->getJD0() + section.times[i]/86400.0;
			std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonBarycenterCS(jd);
			std::array< std::array<double,3>, 3> CS = {frame[1],frame[2],frame[3]};
			std::array<double,3>& origin = frame[0];
			std::array<double,3> pos = {x[0] - origin[0],x[1] - origin[1],x[2] - origin[2]};
			pos = Math::mult(CS,pos);
			fprintf (pFile, format,section.times[i], x[0], x[1], x[2], pos[0], pos[1], pos[2]);
		}
	}
	fclose (pFile);
}

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

void patch_section(Section& section, const std::array<double,6>& xp){

	Matrix<3,4> L;
	Vector<4> u;
	Vector<3> b;
	
	double stepSize = 0.1;
	double norm_old = -1;
	for (int iteration = 0 ; iteration < 20;iteration++) {
		section.compute_states();
		section.compute_STM2();
		
		const Matrix<6,6>& STM = section.STM;
		const std::array<double,6>& xf = section.states.back();
		
		double norm = 0;
		for(uint_fast8_t i = 0; i < 3; i++) {
			for(uint_fast8_t j = 0; j < 3; j++) {
				L[i][j] = STM[i][j+3];
			}
			L[i][3] = xf[i + 3];
			b.data[i] = xp[i] - xf[i];
			norm += b[i]*b[i];
		}
		
		if(norm < 1e-6) {
			break;
		} 

		Matrix<4,3> LT = L.get_transpose();
		Matrix<3,3> A = L*LT;
		try {
			Matrix<3,3>::solve(A,b);
			u = LT*b;
			// u = LT*(MatrixX::invert(A)*b);
		} catch (...) {
			std::cout << "error\n";
			break;
		}

		if(iteration == 0) {
			norm_old = norm;
		}
		
		stepSize *= norm_old/norm;
		norm_old = norm;

		if(stepSize > 1) {
			stepSize = 1;
		} 
		
		if(stepSize < 1e-3) {
			stepSize = 1e-3;
		}

		section.initial_state[3] += stepSize*u[0];
		section.initial_state[4] += stepSize*u[1];
		section.initial_state[5] += stepSize*u[2];
		section.t_final += stepSize*u[3];
	
	}
}


void patch(std::vector<Section>& sections){
	std::cout << "patching...\n";
	const uint_fast16_t nThreads = sections.size() - 1;
	std::vector<std::thread> threads;
	threads.reserve(nThreads);
	
	for(uint_fast16_t i = 0; i < nThreads; i++){
		threads.emplace_back(patch_section, std::ref(sections[i]), std::ref(sections[i+1].initial_state));
	}
	
	for(uint_fast16_t i = 0; i < nThreads; i++){
		threads[i].join();
	}
	
	for (uint_fast16_t section = 1; section <= nThreads; section++) {
		sections[section].t_start = sections[section-1].t_final;
	}
	std::cout << "patched.\n";
}

void minimize_velocity_section(Section& section, const std::array<double,6>& x_front){
	const double deltaX = 0.1;
	const double deltaX_half = deltaX*0.5;
	const double deltaT = 60;
	const double deltaT_half = deltaT*0.5;
	
	double dvx,dvy,dvz,dv1,dv0;
	double dvdx[4];
	std::array<double,6> x_back;
	
	int iter = 0;
	while(iter++ < 20) {
		for(uint_fast8_t i = 0; i < 3; i++){
			section.initial_state[i] -= deltaX_half;
			section.compute_states();
			
			x_back  = section.states.back();
			
			dvx = (x_front[3] - x_back[3]);
			dvy = (x_front[4] - x_back[4]);
			dvz = (x_front[5] - x_back[5]);
			
			dv0 = dvx*dvx + dvy*dvy + dvz*dvz;
			
			section.initial_state[i] += deltaX;
			section.compute_states();
			
			x_back  = section.states.back();
			
			dvx = (x_front[3] - x_back[3]);
			dvy = (x_front[4] - x_back[4]);
			dvz = (x_front[5] - x_back[5]);
			
			dv1 = dvx*dvx + dvy*dvy + dvz*dvz;
			
			section.initial_state[i] -= deltaX_half;
			
			dvdx[i] = (dv1 - dv0)/deltaX;
		}
		
		section.t_start -= deltaT_half;
		section.compute_states();
		
		x_back  = section.states.back();
		
		dvx = (x_front[3] - x_back[3]);
		dvy = (x_front[4] - x_back[4]);
		dvz = (x_front[5] - x_back[5]);
		
		dv0 = dvx*dvx + dvy*dvy + dvz*dvz;
		
		section.t_start += deltaT;
		section.compute_states();
		
		x_back  = section.states.back();
		
		dvx = (x_front[3] - x_back[3]);
		dvy = (x_front[4] - x_back[4]);
		dvz = (x_front[5] - x_back[5]);
		
		dv1 = dvx*dvx + dvy*dvy + dvz*dvz;
		
		section.t_start -= deltaT_half;
		
		dvdx[3] = (dv1 - dv0)/deltaT;
		
		section.compute_states();
		
		dvx = (x_front[3] - x_back[3]);
		dvy = (x_front[4] - x_back[4]);
		dvz = (x_front[5] - x_back[5]);
		
		dv0 = dvx*dvx + dvy*dvy + dvz*dvz;
		
		// we should change initial state as little as possible so
		// dl = dx^2 + dy^2 + dz^2
		// minimize dl
		// subject to: dv + dx*dvdx + dy*dvdy + dz*dvdz = 0
		// 2*dx = gamma*dvdx
		// 2*dy = gamma*dvdy
		// 2*dz = gamma*dvdz
		
		double d = 0.5*dv0/(dvdx[0]*dvdx[0] + dvdx[1]*dvdx[1] + dvdx[2]*dvdx[2]);
		for(uint_fast8_t i = 0; i < 3; i++){
			section.initial_state[i] -= dvdx[i]*d;
		}

		// change dt with regular newton iteration
		double dt = 0.0*dv0/dvdx[3];
		section.t_start -= dt;
		
		section.compute_states();
		
		x_back  = section.states.back();
		
		dvx = (x_front[3] - x_back[3]);
		dvy = (x_front[4] - x_back[4]);
		dvz = (x_front[5] - x_back[5]);
		
		dv1 = dvx*dvx + dvy*dvy + dvz*dvz;
		
		if(dv1 < 1e-6) {
			break;
		}
	}
}

double minimizeDV_Gradient(std::vector<Section>& sections, const EarthMoonSun* dynamics){
	const double deltaX = 1;
	const double deltaX_half = deltaX*0.5;
	const double deltaT = 60;
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
	while(iter++ < 4) {
		for (uint_fast16_t section = 0; section < n1; section++) {
			Section& section_i = sections[section];
			
			const std::array<double,6>& x_front = sections[section+1].initial_state;

			if(section > 0){
				xf_last = sections[section-1].states.back();
			} else {
				xf_last[3] = 0;
				xf_last[4] = 0;
				xf_last[5] = 0;
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
				std::cout << B.to_string() << std::endl;
				Matrix<7,7>::solve(B,x);
				std::cout << B.to_string() << std::endl;
				std::cout << x.to_string() << std::endl;
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

double minimizeDV_Genetic(std::vector<Section>& sections, const EarthMoonSun* dynamics, const double& dx){
	return 0;
}

double minimizeDV(std::vector<Section>& sections, const EarthMoonSun* dynamics){
	const uint_fast16_t nSections = sections.size();
	std::cout << "Minimizing..." << std::endl;
	
	// Recompute STMS, might not be necessary
	for (uint_fast16_t section = 0; section < nSections; section++) {
		sections[section].compute_states();
		sections[section].compute_STM2();
	}
	
	Matrix<3,3> Apf;
	Matrix<3,3> Bpf;
	Matrix<3,3> Cpf;
	Matrix<3,3> Dpf;
	Matrix<3,3> Ap0;
	Matrix<3,3> Bp0;
	Matrix<3,3> Cp0;
	Matrix<3,3> Dp0;
	Matrix<3,3> M0;
	Vector<3> Mt0;
	Matrix<3,3> Mp;
	Vector<3> Mtp;
	Matrix<3,3> Mf;
	Vector<3> Mtf;
	
	Vector<3> vp_back;
	Vector<3> vp_front;
	Vector<3> ap_back;
	Vector<3> ap_front;
	
	MatrixX M(3*nSections-3, 4*nSections + 4,(char)0);
	MatrixX dvp(M.nRows,1);
	MatrixX dr(M.nCols,1);
	MatrixX dv(nSections-1,1);
		
	for (uint_fast16_t section = 0; section < nSections-1; section++) {
		uint_fast16_t colStart = section*4;
		uint_fast16_t rowStart = section*3;
		
		Ap0 = sections[section].STM.slice<3,3>(0,0);
		Bp0 = sections[section].STM.slice<3,3>(0,3);
		Cp0 = sections[section].STM.slice<3,3>(3,0);
		Dp0 = sections[section].STM.slice<3,3>(3,3);
		
		Matrix<6,6> STM_front(sections[section+1].STM);
		Matrix<6,6> STM_reverse = MatrixX::LUPInvert(STM_front);
		
		Apf = STM_reverse.slice<3,3>(0,0);
		Bpf = STM_reverse.slice<3,3>(0,3);
		Cpf = STM_reverse.slice<3,3>(3,0);
		Dpf = STM_reverse.slice<3,3>(3,3);

		const std::array<double,6>& x_back  = dynamics->get_state_rate(sections[section].states.back(), sections[section].times.back());
		const std::array<double,6>& x_front = dynamics->get_state_rate(sections[section+1].initial_state, sections[section+1].t_start);
		
		dvp.data[rowStart] = (x_front[0] - x_back[0]);
		dvp.data[rowStart+1] = (x_front[1] - x_back[1]);
		dvp.data[rowStart+2] = (x_front[2] - x_back[2]);

		for(uint_fast8_t i = 0; i < 3; i++){
			vp_back[i] = x_back[i];
			vp_front[i] = x_front[i];
			ap_back[i] = x_back[i+3];
			ap_front[i] = x_front[i+3];
		}

		const Matrix<3,3>& tmp0 = Dp0*MatrixX::invert(Bp0);
		const Matrix<3,3>& tmpf = Dpf*MatrixX::invert(Bpf);
		
		M0 = tmp0*Ap0 - Cp0;
		Mt0 = ap_back - tmp0*vp_back;

		Mf = Cpf - tmpf*Apf;
		Mtf = tmpf*vp_front - ap_front;
		
		for(uint_fast8_t i = 0; i < 3; i++){
			Mtp.data[i] = - Mt0.data[i] - Mtf.data[i];				
		}
		for(uint_fast8_t i = 0; i < 9; i++){
			Mp.data[i] = tmpf.data[i] - tmp0.data[i];
		}

		for(uint_fast8_t i = 0; i < 3; i++){
			double* row = M[rowStart + i];
			for(uint_fast8_t j = 0; j < 3; j++){
				row[colStart + j] = M0[i][j];
				row[colStart + 4 + j] = Mp[i][j];
				row[colStart + 8 + j] = Mf[i][j];
			}
			row[colStart + 3] = Mt0(i);
			row[colStart + 7] = Mtp(i);
			row[colStart + 11] = Mtf(i);
		}
	}
	
	//std::cout << "DV: " << std::endl << dv.to_string() << std::endl;
	
	MatrixX MT = M.get_transpose();
	MatrixX MMT = M*MT;
	
	Math::LUPSolve(MMT.rows,dvp.data,MMT.nRows);

	for(std::size_t i = 0; i < dr.nRows; i++){
		dr.data[i] = Math::dot(MT[i],dvp.data, MT.nCols);
	}
	
	double norm = 0;
	for (uint_fast16_t section = 0; section < nSections; section++) {
		uint_fast16_t col = section*4;
		for(uint_fast8_t i = 0; i < 3; i++){
			sections[section].initial_state[i] -= dr.data[col + i];
			norm += dr.data[col + i]*dr.data[col + i];
		}
		sections[section].t_start -= dr.data[col + 3];
	}
	// double check times
	for (uint_fast16_t section = 0; section < nSections-1; section++) {
		sections[section].t_final = sections[section+1].t_start;
	}	
	
	return norm;
}

void test(const std::size_t& nSections) {
	
	const double jd0 = Util::getJDFromUTC(2024,4,15,0,0,0);
	std::cout << jd0 << std::endl;	
	EarthMoonSun* dynamics = new EarthMoonSun(jd0);
	
	std::cout << std::setprecision(12);

	// init sections
	std::cout << "Initializing Segments" << std::endl;
	double T = 0;
	const double dT = Section::SECTION_DAYS*86400;
	std::vector<Section> sections(nSections,dynamics);
	for (std::size_t section = 0; section < nSections; section++) { 
	
		double jd = jd0 + T/86400.0;
	
		std::array<double,3> e = dynamics->earth->getPos(jd);
		std::array<double,3> m = dynamics->moon->getPos(jd);
		std::array<double,3> r = {e[0] - m[0],e[1] - m[1],e[2] - m[2]};
		double sma = sqrt(Math::dot(r,r));

		CR3BP cr3bp = CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);
		
		sections[section].initial_state = convert(&cr3bp, dynamics, cr3bp.getHaloInitialState_3rd(10000,0,T,1),jd);
		sections[section].t_start = T;
		T += dT;
		sections[section].t_final = T;
	}
	
	double t_final = T + 3600;
	
	double time = 0;
	std::vector<double> t;
	std::vector<std::array<double,3> > xE;
	std::vector<std::array<double,3> > xM;
	while(time < t_final){
		double jd_t = jd0 + time/86400;
		t.push_back(time);
		xE.push_back(dynamics->earth->getPos(jd_t));
		xM.push_back(dynamics->moon->getPos(jd_t));
		time += 3600;
	}
	
	printOut(t,xE,"test_earth");
	printOut(t,xM,"test_moon");
	
	std::cout << "Printing OG" << std::endl;
	for (Section& section : sections) { 
		section.compute_states();
	}
	printOut(dynamics,sections,"test_orbit");

	int iter = 0;
	double old_error = 1e300;
	double error_min = nSections*1;
	while(iter++ < 5) {

		patch(sections);
		
		try {
			double err = minimizeDV_Gradient(sections,dynamics);
			std::cout << "Mag Error: " << err << std::endl;
			if(fabs(err - old_error)/old_error < 5e-3 || err < error_min){
				//break;
			}
			old_error = err;
		} catch(...) {
			std::cout << "Degenerate Matrix" << std::endl;
			break;
		}

	}
	
	// Final patch
	//patch(sections);
	for (Section& section : sections) { 
		section.compute_states();
	}
	
	std::cout << "printing" << std::endl;
	printOut(dynamics,sections,"test_orbit2");
	
	double deltaV = 0;
	for (uint_fast16_t section = 0; section < nSections-1; section++) {

		const std::array<double,6>& x_back  = sections[section].states.back();
		const std::array<double,6>& x_front = sections[section+1].initial_state;
		
		double dvx = (x_front[3] - x_back[3]);
		double dvy = (x_front[4] - x_back[4]);
		double dvz = (x_front[5] - x_back[5]);
		deltaV += sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
	}
	std::cout << "approx dv: " << deltaV << std::endl;
}

void test2(){
	
	CR3BP cr3bp = CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,385000);
	
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(&cr3bp);
	std::array<double,6> x = cr3bp.getHaloInitialState_3rd(10000,0,120*86400,1);
	ode.recording.set_record_interval(0.01);
	ode.set_timestep(0.0001);
	
	ode.run(x,4);
	
	printOut(ode.recording,"test_orbit3");
	std::cout << "done";
}

int main(int argc, char* argv[]) {

	int nSections = 30;
	if(argc > 1) {
		nSections = (int)(std::stod(argv[1])/Section::SECTION_DAYS);
	}
	test(nSections);
	
	return 0;
}
