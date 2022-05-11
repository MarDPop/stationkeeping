#include "Dynamics.h"
#include "ODE.h"
#include "Util.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define MAX_ITERATIONS 10

void printOut(const std::array<double,6>& corrected){
	FILE * pFile;

	pFile = fopen ("initial_state","w");
	for (int i = 0 ; i < 6; i++) {
		fprintf (pFile, "%16.14f ",corrected[i]);
	}
	fclose (pFile);
}

void printOut(const Recording<6>&  record){
	FILE * pFile;

	pFile = fopen ("orbit","w");
	const int n = record.number_entries();
	std::string formatStr = "%12.6f %16.14f %16.14f %16.14f\n";
	const char* format = formatStr.c_str();
	for (int i = 0 ; i < n; i++) {
		const std::array<double,6>& state = record.get(i);
		fprintf (pFile, format,record.time_at(i), state[0], state[1], state[2]);
	}
	fclose (pFile);
}

void runFull(EarthMoonSun* dynamics, const double& Az, const double& jd0, const std::array<double,6>& initial_state, double tol, double recording_interval){
	
	dynamics->setJD(jd0,jd0+150,0.5);
	
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(dynamics);
	ode.recording.set_record_interval(recording_interval);
	ode.set_timestep(recording_interval*1e-3);
	ode.stop = [](const std::array<double,6>& x,const double& t){
		return x[1] < 0;
	};
	
	std::array<double,6> x = initial_state;
	
	std::cout << "Running" <<std::endl;
	double** A = Math::zeros(6);
	double** STM = Math::zeros(6);
	double** dSTM = Math::zeros(6);
	for(int iter = 0; iter < MAX_ITERATIONS;iter++) {
		ode.recording.clear();
		ode.run(x,100);
		
		std::array<double,6> xf = ode.get_state();
		
		std::array<double,6> dxf = dynamics->get_state_rate(xf,ode.get_time());
		
		double dtf = -xf[1]/dxf[1];
		double tf = ode.get_time() + dtf;
		
		for(int i = 0; i < 6; i++){
			xf[i] += dxf[i]*dtf;
		}
		
		std::cout << "half period z velocity (m/s): " << xf[5] << std::endl;
		
		double e = fabs(xf[5]);
        
        if (e < tol)
            break;
		
		double jdf = jd0 + tf/86400;
		std::array<double,3> earth = dynamics->getEarth_EMB(jdf);
		std::array<double,3> moon = dynamics->getMoon_EMB(jdf);
		std::array<double,3> d = {earth[0] - moon[0],earth[1] - moon[1],earth[2] - moon[2]};
		double sma = sqrt(Math::dot(d,d));

		CR3BP cr3bp = CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);
		
		cr3bp.getInitialState(Az,0,tf,false);
    
        //Integrate State Transition Matrix
        Math::identity(STM,6);
        int n = ode.recording.number_entries()-1;
        for(int i = 0; i < n;i++) {
			const std::array<double,6>& xi = ode.recording.get(i);
			double dt = ode.recording.time_at(i+1) - ode.recording.time_at(i);
			double jd = jd0 + ode.recording.time_at(i)/86400;
            dynamics->getA(xi,jd,A);
            Math::mult(A,STM,dSTM,6);
			Math::mult(dSTM,dt,6);
            Math::add(STM,dSTM,6);
		}	
		dynamics->getA(ode.recording.get(n),jd0 + ode.recording.time_at(n)/86400,A);
		Math::mult(A,STM,dSTM,6);
		Math::mult(dSTM,dtf,6);
		Math::add(STM,dSTM,6);	
    
        double stepSize = 1/(1+e*5);
    
		double dx0 = 0;
		double dz0 = 0;
		double dyd0 = 0;
		x[0] += stepSize*dx0;
		x[2] += stepSize*dz0;
        x[4] += stepSize*dyd0;
    
        std::cout << "corrected initial state: " << std::endl;
		
		for(int i = 0; i < 6;i++){
			std::cout << x[i] << " ";
		}
		std::cout << std::endl;
	}
	
	Math::del(A,6);
	Math::del(STM,6);
	Math::del(dSTM,6);
	
	printOut(x);
	printOut(ode.recording);
}

int main(int argc, char* argv[]) { 
	
	double Az = 10000;
	double tol = 1e-6;
	double recording_interval = 3600;
	
	std::array<double,6> initial_state;
	double jd0 = 0; 
	
	bool useInitialState = false;
	bool estimateInitialState = false;
	for(int arg = 1; arg < argc;arg++){
		if(argv[arg][0] == '-'){
			switch(argv[arg][1]){
				case 'j':
					if (arg + 1 < argc){
						jd0 = std::stod(argv[arg+1]);
					} else {
						std::cout << "julian date must be specificed";
						return 1;
					}
					break;
				case 'g':
					if (arg + 6 < argc){
						int year = std::stoi(argv[arg+1]);
						int mon = std::stoi(argv[arg+2]);
						int day = std::stoi(argv[arg+3]);
						int hour = std::stoi(argv[arg+4]);
						int min = std::stoi(argv[arg+5]);
						int sec = std::stoi(argv[arg+6]);
						jd0 = Util::getJDFromUTC(year,mon,day,hour,min,sec);
					} else {
						std::cout << "julian date must be specificed";
						return 1;
					}
					break;
				case 'o':
					if (arg + 1 < argc){
						recording_interval = std::stod(argv[arg+1]);
					} 
					break;
				case 'i':
					useInitialState = true;
					if (arg + 6 < argc){
						for(int i = 0; i < 6; i++)
							initial_state[i] = std::stod(argv[arg+i+1]);
					} else {
						std::cout << "Need 6 state positions";
						return 1;
					}
					break;
				case 'e':
					estimateInitialState = true;
					if (arg + 1 < argc){
						Az = std::stod(argv[arg+1]);
					} else {
						std::cout << "Need Halo Size in km";
						return 1;
					}
					break;
				case 't':
					if (arg + 1 < argc){
						tol = std::stod(argv[arg+1]);
					} else {
						std::cout << "Specify tolerance in km/s";
						return 1;
					}
					break;
				case 's':
					if (arg + 1 < argc){
						tol = std::stod(argv[arg+1]);
					} else {
						std::cout << "Specify tolerance in km/s";
						return 1;
					}
					break;
				case 'h':
					std::cout << "-o <Optional/interval>\n";
					std::cout << "-i <x> <y> <z> <u> <v> <w> initial state\n";
					std::cout << "-s Use solar exclusion\n";
					std::cout << "-e <Az> estimate state from halo\n";
					std::cout << "-t <Xtol> <Vtol> tolerance for \n";
					break;
			}
		}
	}
	
	if(!estimateInitialState && !useInitialState) {
			std::cout << "need to define intial state";
			return 2;
		}
	
	if(jd0 < 2400000) {
		std::cout << "need to define julian date";
		return 3;
	}
	
	EarthMoonSun* dynamics = new EarthMoonSun();
	
	if(estimateInitialState) {		
		std::array<double,3> e = dynamics->getEarth_EMB(jd0);
		std::array<double,3> m = dynamics->getMoon_EMB(jd0);
		
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
		std::cout << "jd: " << jd0 << std::endl;
		std::cout << "Earth-Moon Distance (km): " << sma << std::endl;
		
		CR3BP cr3bp = CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);
		
		std::array<double,6> state_guess = cr3bp.getInitialState(Az,0,0,false);
		
		std::cout << "State Normalized CR3BP Guess" << std::endl;
		for(std::size_t i = 0 ; i < 6;i++){
			std::cout << state_guess[i] << " ";
		}
		std::cout << std::endl;
		
		std::array<double,6> x = cr3bp.convert_state_to_inertial(state_guess);

		initial_state = dynamics->cr3bp_to_embj2000(x,jd0);
	}

	runFull(dynamics,Az,jd0,initial_state,tol,recording_interval);

	delete dynamics;
	return 0;
}
