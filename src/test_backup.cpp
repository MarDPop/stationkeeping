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
#include "../Eigen/Dense"

#define MAX_ITERATIONS 10

//typedef Eigen::Matrix<double,6,6> Matrix6d;
//typedef Eigen::Matrix<double,6,1> Vector6d;

struct Section {
	std::array<double,6> initial_state;
	double t_start, t_final;
	std::vector<double> times;
	std::vector< std::array<double,6> > states;
	Matrix<6,6> STM;
	
	void compute_states(ODE_RK4<6>* ode){
		ode->recording.clear();
		ode->set_time(this->t_start);
		ode->run(this->initial_state,this->t_final);
		this->times = ode->recording.get_times();
		this->states = ode->recording.get_states();
		this->times.push_back(ode->get_time());
		this->states.push_back(ode->get_state());
	}
	
	void compute_STM(Dynamics<6>* dynamics){
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
			a = dynamics->get_state_rate(x, this->times[i]);
			x[0] += 1e-3;
			ax = dynamics->get_state_rate(x, this->times[i]);
			x[0] -= 1e-3;
			x[1] += 1e-3;
			ay = dynamics->get_state_rate(x, this->times[i]);
			x[1] -= 1e-3;
			x[2] += 1e-3;
			az = dynamics->get_state_rate(x, this->times[i]);
			
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
	
	void compute_STM2(EarthMoonSun* dynamics){
		Matrix<6,6> A;
		this->STM.set_identity();
		const int n = this->times.size()-1;
        for(int i = 0; i < n;i++) {
			double dt = this->times[i+1] - this->times[i];
			double jd = dynamics->getJD0() + this->times[i]/86400;
            dynamics->getA(this->states[i],jd,A);
			STM += A*STM*dt;
		}
	}
	
	void compute_STM3(EarthMoonSun* dynamics){
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
			double jd = dynamics->getJD0() + this->times[i]/86400;
			Gs.push_back(std::move(dynamics->getG(this->states[i],jd)));
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

void printOut(const Recording<6>&  record){
	FILE * pFile;

	pFile = fopen ("test_orbit","w");
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
	
	for(Section section : sections) {
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

void testKepler(){
	double deg2rad = 0.017453292519943295769236907684;
	std::array<double,6> oe = {4.691931220181541E+03,3.752109909226587E-02,5.244919734790719*deg2rad,7.962340298182825E+01*deg2rad,3.163784072969009E+02*deg2rad,266.55391327771815*deg2rad };
	std::array<double,3> pos = OrbitalElements::kepler2position(oe);
	std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	double TA = OrbitalElements::trueAnomalyFromMeanAnomaly(2.708520044993519E+02*deg2rad,3.752109909226587E-02);
	double TA2 = OrbitalElements::trueAnomalyFromMeanAnomalyApprox(2.708520044993519E+02*deg2rad,3.752109909226587E-02);
	std::cout << (TA/deg2rad + 360) << " VS " << (TA2/deg2rad + 360) << std::endl;
}

std::array<double,6> convert(CR3BP* cr3bp, EarthMoonSun* dynamics, const std::array<double,6>& state_guess, const double& jd){
	std::array<double,6> x = cr3bp->convert_state_to_inertial(state_guess);
	
	std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonBarycenterCS(jd);
	
	std::array< std::array<double,3>, 3> CS = {frame[1],frame[2],frame[3]};
	std::array< std::array<double,3>, 3> CST = Math::transpose(CS);
	std::array<double,3>& origin = frame[0];
	std::array<double,3> pos = {origin[0] + x[0],origin[1] + x[1],origin[2] + x[2]};
	std::array<double,3> vel = {x[3],x[4],x[5]};
	
	pos = Math::mult(CST,pos);
	vel = Math::mult(CST,vel);

	x[0] = pos[0];
	x[1] = pos[1];
	x[2] = pos[2];
	x[3] = vel[0];
	x[4] = vel[1];
	x[5] = vel[2];
	
	return x;
}

void patch(std::vector<Section>& sections, ODE_RK4<6>* ode){
	// https://webthesis.biblio.polito.it/6898/1/tesi.pdf
	// https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1988CeMec..41..107H&defaultprint=YES&filetype=.pdf
	const std::size_t nSections = sections.size();
	std::cout << "Patching..." << std::endl;
	
	EarthMoonSun* dynamics = (EarthMoonSun*)ode->get_dynamics();
	
	std::vector<bool> finished(nSections,false);
	std::vector<double> stepSize(nSections,0.5);
	std::vector<double> norms(nSections,400.0);
	
	Matrix<3,4> L;
	Vector<4> u;
	Vector<3> b;
		
	for (int iteration = 0 ; iteration < 20;iteration++) {
		std::cout << "iteration: " << iteration << std::endl;
		int completed = 0;
		for (std::size_t section = 0; section < nSections-1; section++) {
			
			if (finished[section]) {
				completed++;
				continue;
			}
			sections[section].compute_states(ode);
			sections[section].compute_STM2(dynamics);
			
			const Matrix<6,6>& STM = sections[section].STM;
			const std::array<double,6>& xf = sections[section].states.back();
			const std::array<double,6>& xp = sections[section+1].initial_state;
			
			double norm = 0;
			for(std::size_t i = 0; i < 3; i++) {
				for(std::size_t j = 0; j < 3; j++) {
					L[i][j] = STM[i][j+3];
				}
				L[i][3] = xf[i + 3];
				b[i] = xp[i] - xf[i];
				norm += b[i]*b[i];
			}
			
			if(norm < 0.0001) {
				finished[section] = true;
				continue;
			} 
			
			double dnorm = (norm - norms[section])/stepSize[section];
			norms[section] = norm;
			stepSize[section] -= 0.9*norm/dnorm;

			if(stepSize[section] > 1) {
				stepSize[section] = 1;
			} 
			
			if(stepSize[section] < 1e-4) {
				stepSize[section] = 1e-4;
			}
			
			Matrix<4,3> LT = L.get_transpose();
			Matrix<3,3> A = L*LT;
			Matrix<3,3>::solve(A,b);
			u = LT*b;

			sections[section].initial_state[3] += stepSize[section]*u[0];
			sections[section].initial_state[4] += stepSize[section]*u[1];
			sections[section].initial_state[5] += stepSize[section]*u[2];
			sections[section].t_final += stepSize[section]*u[3];
		}
		
		std::cout << completed << " Sections Complete" << std::endl;
		
		if(completed == nSections-1) {
			break;
		}
	}
	
}

double minimizeDV(std::vector<Section>& sections, ODE_RK4<6>* ode){
	const std::size_t nSections = sections.size();
	std::cout << "Minimizing..." << std::endl;
	
	EarthMoonSun* dynamics = (EarthMoonSun*)ode->get_dynamics();
	
	// Recompute STMS, might not be necessary
	for (std::size_t section = 0; section < nSections; section++) {
		sections[section].compute_states(ode);
		sections[section].compute_STM2(dynamics);
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
		
	for (std::size_t section = 0; section < nSections-1; section++) {
		std::size_t colStart = section*4;
		std::size_t rowStart = section*3;
		
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
		const std::array<double,6>& x_front = dynamics->get_state_rate(sections[section+1].states[0], sections[section+1].times[0]);
		
		dvp.data[rowStart] = (x_front[0] - x_back[0]);
		dvp.data[rowStart+1] = (x_front[1] - x_back[1]);
		dvp.data[rowStart+2] = (x_front[2] - x_back[2]);
		dv.data[section] = sqrt(dvp.data[rowStart]*dvp.data[rowStart] + dvp.data[rowStart+1]*dvp.data[rowStart+1] + dvp.data[rowStart+2]*dvp.data[rowStart+2]);
		for(std::size_t i = 0; i < 3; i++){
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
		
		for(std::size_t i = 0; i < 3; i++){
			Mtp.data[i] = - Mt0.data[i] - Mtf.data[i];				
		}
		for(std::size_t i = 0; i < 9; i++){
			Mp.data[i] = tmpf.data[i] - tmp0.data[i];
		}

		for(std::size_t i = 0; i < 3; i++){
			double* row = M[rowStart + i];
			for(std::size_t j = 0; j < 3; j++){
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
		dr.data[i] = 0.5*Math::dot(MT[i],dvp.data, MT.nCols);
	}
	
	double norm = 0;
	for (std::size_t section = 0; section < nSections; section++) {
		std::size_t col = section*4;
		for(std::size_t i = 0; i < 3; i++){
			sections[section].initial_state[i] -= dr.data[col + i];
			norm += dr.data[col + i]*dr.data[col + i];
		}
		sections[section].t_start -= dr.data[col + 3];
	}
	// double check times
	for (std::size_t section = 0; section < nSections-1; section++) {
		sections[section].t_final = sections[section+1].t_start;
	}
	
	std::cout << "Mag Error: " << norm << std::endl;	
	
	return norm;
}

void test3(const std::size_t& nSections) {
	
	EarthMoonSun* dynamics = new EarthMoonSun();
	double jd0 = Util::getJDFromUTC(2021,1,1,0,0,0);
	dynamics->setJD(jd0,jd0 + 60,0.1);
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(dynamics);
	ode.recording.set_record_interval(300);
	ode.set_timestep(5);
	
	std::cout << std::setprecision(12);
	
	// init sections
	double day = 0;
	const double deltaDay = .75;
	std::vector<Section> sections(nSections);
	for (std::size_t section = 0; section < nSections; section++) { 
	
		double jd = jd0 + day;
	
		std::array<double,3> e = dynamics->getEarth_EMB(jd);
		std::array<double,3> m = dynamics->getMoon_EMB(jd);
		std::array<double,3> r = {e[0] - m[0],e[1] - m[1],e[2] - m[2]};
		double sma = sqrt(Math::dot(r,r));

		CR3BP cr3bp = CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);
		
		sections[section].initial_state = convert(&cr3bp, dynamics, cr3bp.getHaloInitialState_3rd(10000,0,day*86400,1),jd);
		sections[section].t_start = day*86400;
		sections[section].t_final = (day+deltaDay)*86400;
		
		ode.set_time(sections[section].t_start);
		ode.run(sections[section].initial_state,sections[section].t_final);
		
		sections[section].times = ode.recording.get_times();
		sections[section].states = ode.recording.get_states();
		sections[section].times.push_back(ode.get_time());
		sections[section].states.push_back(ode.get_state());

		ode.recording.clear();
		day += deltaDay;
	}
	
	double t_final = day*86400 + 3600;
	
	double time = 0;
	std::vector<double> t;
	std::vector<std::array<double,3> > xE;
	std::vector<std::array<double,3> > xM;
	while(time < t_final){
		double jd_t = jd0 + time/86400;
		t.push_back(time);
		xE.push_back(dynamics->getEarth_EMB(jd_t));
		xM.push_back(dynamics->getMoon_EMB(jd_t));
		time += 900;
	}
	
	printOut(t,xE,"test_earth");
	printOut(t,xM,"test_moon");
	
	std::cout << "Printing OG" << std::endl;
	printOut(dynamics,sections,"test_orbit");
	
	int iter = 0;
	const double tol = nSections*1;
	while(iter++ < 16) {

		patch(sections,&ode);
		
		if(minimizeDV(sections,&ode) < tol){
			break;
		}
	}
	
	// Final patch
	patch(sections,&ode);
	
	std::cout << "printing" << std::endl;
	printOut(dynamics,sections,"test_orbit2");

}

void test2(){
	
	EarthMoonSun* dynamics = new EarthMoonSun();
	double jd = Util::getJDFromUTC(2021,1,1,0,0,0);
	dynamics->setJD(jd,jd+150,0.5);
	
	std::array<double,3> e = dynamics->getEarth_EMB(jd);
	std::array<double,3> m = dynamics->getMoon_EMB(jd);
	
	std::cout << "Earth Position (km)" << std::endl;
	for(std::size_t i = 0 ; i < 3;i++){
		std::cout << e[i] << " ";
	}
	std::cout << std::endl;
	
	std::cout << "Moon Position (km)" << std::endl;
	for(std::size_t i = 0 ; i < 3;i++){
		std::cout << m[i] << " ";
	}
	std::cout << std::endl;
	
	std::array<double,3> r;
	for(std::size_t i = 0; i < 3;i++) {
		r[i] = e[i] - m[i];
	}
	double sma = sqrt(Math::dot(r,r));
	std::cout << std::setprecision(15);
	std::cout << "jd: " << jd << std::endl;
	std::cout << "Earth-Moon Distance (km): " << sma << std::endl;
	
	CR3BP* cr3bp = new CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);
	
	std::array<double,6> state_guess = cr3bp->getHaloInitialState_3rd(10000,0,0,1);
	
	std::cout << "State Normalized CR3BP Guess" << std::endl;
	for(std::size_t i = 0 ; i < 6;i++){
		std::cout << state_guess[i] << " ";
	}
	std::cout << std::endl;
	
	std::array<double,6> x = cr3bp->convert_state_to_inertial(state_guess);
	
	std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonBarycenterCS(jd);
	std::cout << "frame" << std::endl;
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 3;j++){
			std::cout<< frame[i][j] << " ";
		}
		std::cout << std::endl;
	}
	
	std::cout << "State CR3BP Guess" << std::endl;
	for(std::size_t i = 0 ; i < 6;i++){
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;
	
	std::array< std::array<double,3>, 3> CS = {frame[1],frame[2],frame[3]};
	std::array< std::array<double,3>, 3> CST = Math::transpose(CS);
	std::array<double,3>& origin = frame[0];
	std::array<double,3> pos = {origin[0] + x[0],origin[1] + x[1],origin[2] + x[2]};
	std::array<double,3> vel = {x[3],x[4],x[5]};
	
	pos = Math::mult(CST,pos);
	vel = Math::mult(CST,vel);

	x[0] = pos[0];
	x[1] = pos[1];
	x[2] = pos[2];
	x[3] = vel[0];
	x[4] = vel[1];
	x[5] = vel[2];
	
	std::cout << "State Guess Inertial (km km/s)" << std::endl;
	for(std::size_t i = 0 ; i < 6;i++){
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;
	
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(dynamics);
	ode.recording.set_record_interval(600);
	ode.set_timestep(1);
	
	double deltaT = 86400*2;

	ode.run(x,deltaT);
	
	printOut(ode.recording);
	
	double time = 0;
	std::vector<double> t;
	std::vector<std::array<double,3> > xE;
	std::vector<std::array<double,3> > xM;
	while(time < ode.get_time() + 3600){
		double jd_t = jd + time/86400;
		t.push_back(time);
		xE.push_back(dynamics->getEarth_EMB(jd_t));
		xM.push_back(dynamics->getMoon_EMB(jd_t));
		time += 900;
	}
	
	printOut(t,xE,"test_earth");
	printOut(t,xM,"test_moon");
	
}

void test1(){
	
	//double beta = 3.04036e-6;
	//double a = 1.49598e8;
	//CR3BP* cr3bp = new CR3BP(OrbitalElements::SUN_MU,OrbitalElements::SUN_MU*beta/(1-beta),a);
	//std::array<double,6> state_guess = cr3bp->getHaloInitialState_3rd(125000,0,0,1);
	CR3BP* cr3bp = new CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,385000);
	std::array<double,6> state_guess = cr3bp->getHaloInitialState_3rd(10000,0,86400*1,1);
	
	std::cout << "State Normalized CR3BP Guess" << std::endl;
	for(std::size_t i = 0 ; i < 6;i++){
		std::cout << state_guess[i] << " ";
	}
	std::cout << std::endl;
	
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(cr3bp);
	ode.recording.set_record_interval(0.001);
	ode.set_timestep(1e-6);

	ode.run(state_guess,3);
	
	printOut(ode.recording);	
}

void testInverseMatrix() {
	Matrix<4,4> A; 
	Eigen::Matrix<double,4,4> _A;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			double val = rand() % 100;
			A[i][j] = val;
			_A(i,j) = val;
		}
	}
	Matrix<4,4> Ainv = MatrixX::LUPInvert<4>(A);
	Eigen::Matrix<double,4,4> _Ainv = _A.inverse();
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			std::cout << Ainv[i][j] << " ";
		}
		std::cout << std::endl;
	}
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			std::cout << _Ainv(i,j) << " ";
		}
		std::cout << std::endl;
	}
}

int main(int argc, char* argv[]) {
	//testInverseMatrix();	
	test3(30);
	
	return 0;
}
