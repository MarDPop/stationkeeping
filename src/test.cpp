#include "../include/OrbitalElements.h"
#include "../include/OrbitalDynamics.h"
#include "../include/ODE.h"
#include "../include/Math.h"
#include "../include/Util.h"
#include "../include/Matrix.h"
#include "../include/Simulation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

//#include "../Eigen/Dense"

//typedef Eigen::Matrix<double,6,6> Matrix6d;
//typedef Eigen::Matrix<double,6,1> Vector6d;

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

	std::cout << "Printing OG" << std::endl;
	
	printOut(t,xE,"test_earth");
	printOut(t,xM,"test_moon");

	for (Section& section : sections) { 
		section.compute_states();
	}
	printOut(dynamics,sections,"test_orbit");

	double dvEnd = nSections*0.01;
	const int MAX_ITER = 30;
	for(int iter = 0; iter < MAX_ITER; iter++) {

		Section::minimizeDX(sections);
		
		Section::minimizeDV(sections);
		
		double dv = Section::calcDV(sections);
		std::cout << "DV: " << dv << std::endl;
		if(dv < dvEnd){
			break;
		}

	}
	
	std::cout << "printing" << std::endl;
	printOut(dynamics,sections,"test_orbit2");
}

int main(int argc, char* argv[]) {

	int nSections = 30;
	if(argc > 1) {
		nSections = (int)(std::stod(argv[1])/Section::SECTION_DAYS);
	}
	test(nSections);
	
	return 0;
}
