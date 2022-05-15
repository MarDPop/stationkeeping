#include "../include/OrbitalDynamics.h"
#include "../include/OrbitalElements.h"
#include "../include/ODE.h"
#include "../include/Util.h"
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

//https://www.sciencedirect.com/sdfe/reader/pii/S1000936118300980/pdf
//https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?2008ASPC..394..734S&defaultprint=YES&filetype=.pdf
//https://drum.lib.umd.edu/bitstream/handle/1903/2120/umi-umd-2090.pdf?sequence=1&isAllowed=y
void runCR3BP(const double& MU1, const double& MU2, const double& SMA, const std::array<double,6>& initial_state, double tol, double recording_interval){
	CR3BP* dynamics = new CR3BP(MU1,MU2,SMA);
	
	double mean_motion = sqrt((MU1 + MU2)/(SMA*SMA*SMA));
    double f = SMA*mean_motion;
    tol /= f;
	recording_interval *= mean_motion; 
	
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(dynamics);
	ode.recording.set_record_interval(recording_interval);
	ode.set_timestep(1e-6);
	ode.stop = [](const std::array<double,6>& x,const double& t){
		return x[1] < 0;
	};
	
	std::array<double,6> x = initial_state;
	
	std::cout << "Running" <<std::endl;
	double A00,A11,A10,A01,dx0,dz0,dyd0;
	double** A = Math::zeros(6);
	double** STM = Math::zeros(6);
	double** dSTM = Math::zeros(6);
	bool alternate = true;
	for(int iter = 0; iter < MAX_ITERATIONS;iter++) {
		ode.recording.clear();
		ode.run(x,100);
		
		std::array<double,6> xf = ode.get_state();
		
		std::array<double,6> dxf = dynamics->get_state_rate(xf,ode.get_time());
		
		double dtf = -xf[1]/dxf[1];
		
		double x_dot_f = xf[3] + dxf[3]*dtf;
		double z_dot_f = xf[5] + dxf[5]*dtf;
		
		std::cout << "half period x velocity (m/s): " << x_dot_f*f*1e3 << std::endl;
		std::cout << "half period z velocity (m/s): " << z_dot_f*f*1e3 << std::endl;
		
		double e = (fabs(x_dot_f) + fabs(z_dot_f));
        
        if (e < tol)
            break;
    
        //Integrate State Transition Matrix
        Math::identity(STM,6);
        int n = ode.recording.number_entries()-1;
        for(int i = 0; i < n;i++) {
			const std::array<double,6>& xi = ode.recording.get(i);
			double dt = ode.recording.time_at(i+1) - ode.recording.time_at(i);
            dynamics->getA(xi,A);
            Math::mult(A,STM,dSTM,6);
			Math::mult(dSTM,dt,6);
            Math::add(STM,dSTM,6);
		}	
		dynamics->getA(ode.recording.get(n),A);
		Math::mult(A,STM,dSTM,6);
		Math::mult(dSTM,dtf,6);
		Math::add(STM,dSTM,6);	
    
        //Solve for delta intial state
		
        double xdd = dxf[3]/xf[4];
        double zdd = dxf[5]/xf[4];
    
        double stepSize = 1/(1+e*5);
    
        //alternate between dx0 and dx0 matrices
        if (alternate) {
            A00 = STM[3][0] - xdd*STM[1][0];
            A01 = STM[3][4] - xdd*STM[1][4];
            A10 = STM[5][0] - zdd*STM[1][0];
            A11 = STM[5][4] - zdd*STM[1][4];
            double det = A00*A11 - A10*A01;
            dx0 = (A01*z_dot_f - A11*x_dot_f)/det;
            dyd0 = (A10*x_dot_f - A00*z_dot_f)/det;
            x[0] += stepSize*dx0;
        } else{
            A00 = STM[3][2] - xdd*STM[1][2];
            A01 = STM[3][4] - xdd*STM[1][4];
            A10 = STM[5][2] - zdd*STM[1][2];
            A11 = STM[5][4] - zdd*STM[1][4];
            double det = A00*A11 - A10*A01;
            dz0 = (A01*z_dot_f - A11*x_dot_f)/det;
            dyd0 = (A10*x_dot_f - A00*z_dot_f)/det;
            x[2] += stepSize*dz0;
		}
        x[4] += stepSize*dyd0;
    
        std::cout << "corrected initial state: " << std::endl;
		
		for(int i = 0; i < 6;i++){
			std::cout << x[i] << " ";
		}
		std::cout << std::endl;
    
        alternate = !alternate;
	}
	
	Math::del(A,6);
	Math::del(STM,6);
	Math::del(dSTM,6);
	
	printOut(x);
	printOut(ode.recording);
}

int main(int argc, char* argv[]) { 

	double MU1 = OrbitalElements::EARTH_MU;
	double MU2 = OrbitalElements::MOON_MU;
	double SMA = 385000;
	
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
				case 'o':
					if (arg + 1 < argc){
						recording_interval = std::stod(argv[arg+1]);
					} 
					break;
				case 'm':
					if (arg + 2 < argc){
						MU1 = std::stod(argv[arg+1]);
						MU2 = std::stod(argv[arg+2]);
					} else {
						std::cout << "Need 2 Body Masses";
						return 1;
					}
					break;
				case 'r':
					if (arg + 1 < argc){
						SMA = std::stod(argv[arg+1]);
					} else {
						std::cout << "Specify Distance";
						return 1;
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
				case 'h':
					std::cout << "-f <file>\n";
					std::cout << "-o <Optional/interval>\n";
					std::cout << "-m <mu1> <mu2>\n";
					std::cout << "-r <distance>\n";
					std::cout << "-i <x> <y> <z> <u> <v> <w>\n";
					std::cout << "-e <Az>\n";
					std::cout << "-t <tol>\n";
					break;
			}
		}
	}
	
	if(!estimateInitialState && !useInitialState) {
			std::cout << "need to define intial state";
			return 2;
		}
	
	if(estimateInitialState) {
		CR3BP* dynamics = new CR3BP(MU1,MU2,SMA);
		initial_state = dynamics->getInitialState(Az);
		std::cout << "Initial State\n";
		for(int i = 0; i < 6;i++){
			std::cout << initial_state[i] << " ";
		}
		std::cout << std::endl;
	}
	
	runCR3BP(MU1,MU2,SMA,initial_state,tol,recording_interval);

	return 0;
}
