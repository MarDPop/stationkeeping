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

#define MAX_ITERATIONS 10

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
		this->ode.set_timestep(5);
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

const double Section::SECTION_DAYS = 0.5;

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
	
	const double xL1 = cr3bp->getL1();
	
	std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonBarycenterCS(jd);
	
	std::array<double,3> moon = dynamics->moon->getPos(jd - 0.0005);
	std::array<double,3> earth = dynamics->earth->getPos(jd - 0.0005);
	std::array<double,3> r0 = {moon[0] - earth[0],moon[1] - earth[1],moon[2] - earth[2]};
	
	moon = dynamics->moon->getPos(jd + 0.0005);
	earth = dynamics->earth->getPos(jd + 0.0005);
	std::array<double,3> r1 = {moon[0] - earth[0],moon[1] - earth[1],moon[2] - earth[2]};
	
	double s0 = 0;
	double s1 = 0;
	for(int i = 0; i < 3;i++){
		s0 += r0[i]*r0[i];
		s1 += r1[i]*r1[i];
	}
	double dr = sqrt(s1) - sqrt(s0);
	double dt = 0.001*86400;
	double dL1 = xL1*dr/dt;
	
	//std::cout << jd << ": " << dL1 << std::endl;
	
	std::array< std::array<double,3>, 3> CS = {frame[1],frame[2],frame[3]};
	std::array< std::array<double,3>, 3> CST = Math::transpose(CS);
	
	std::array<double,3> pos = {x[0],x[1],x[2]};
	std::array<double,3> vel = {x[3] + dL1,x[4],x[5]};
	
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
	std::cout << "patched.\n";
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
		dr.data[i] = 0.2*Math::dot(MT[i],dvp.data, MT.nCols);
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
	
	std::cout << "Mag Error: " << norm << std::endl;	
	
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
	double error_min = nSections*100;
	while(iter++ < 5) {

		patch(sections);
		
		for (uint_fast16_t section = 1; section < nSections; section++) {
			sections[section].t_start = sections[section-1].t_final;
		}
		
		try {
			double err = minimizeDV(sections,dynamics);
			if(fabs(err - old_error)/old_error < 5e-3 || err < error_min){
				break;
			}
			old_error = err;
		} catch(...) {
			std::cout << "Degenerate Matrix" << std::endl;
			break;
		}

	}
	
	// Final patch
	patch(sections);
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
