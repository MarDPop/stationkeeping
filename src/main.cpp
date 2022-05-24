#include "../include/Simulation.h"
#include "../include/OrbitalDynamics.h"
#include "../include/OrbitalElements.h"
#include "../include/Util.h"
#include "../include/Config.h"
#include "../include/Math.h"
#include <cmath>
#include <array>

std::array<double,3> convert(CR3BP* cr3bp, EarthMoonSun* dynamics, const std::array<double,3>& pos, const double& jd){
	std::array<double,3> inertial = {pos[0]*cr3bp->sma,pos[1]*cr3bp->sma,pos[2]*cr3bp->sma};
	
	std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonBarycenterCS(jd);

	std::array< std::array<double,3>, 3> CST = {{{frame[1][0],frame[2][0],frame[3][0]},{frame[1][1],frame[2][1],frame[3][1]},{frame[1][2],frame[2][2],frame[3][2]}}};

	return Math::mult(CST,pos);
}

void test(const int& nSections, const double& jd, const double& Az) {
	VectorTable earth = VectorTable("resources/EARTH_EMB_VECTOR.dat");
    VectorTable moon = VectorTable("resources/MOON_EMB_VECTOR.dat");

	std::array<double,3> earth_pos = earth.getPos(jd);
	std::array<double,3> moon_pos = moon.getPos(jd);
	double dx = earth_pos[0] - moon_pos[0];
	double dy = earth_pos[1] - moon_pos[1];
	double dz = earth_pos[2] - moon_pos[2];
	double sma = sqrt(dx*dx + dy*dy + dz*dz);

	CR3BP* cr3bp = new CR3BP(OrbitalElements::EARTH_MU, OrbitalElements::MOON_MU,sma);

	std::array<double,6> x0 = cr3bp->get_halo_initial_state_3rd_order(Az,0,0,1);
	std::array<double,6> xf;
	double tf;
	OrbitComputation::differential_correct_cr3bp(cr3bp,x0,xf,tf,true,1e-10);

	double period = tf*2/cr3bp->mean_motion;
	std::cout << "Period Guess: " << period/Util::JULIAN_DAY << "days vs. " << cr3bp->get_period(Az)/86400 << "days." << std::endl;

	Recording<6> orbit = OrbitComputation::get_cr3bp_halo_orbit(cr3bp,x0,60,period);

	Util::printOut(orbit,"output/CR3BP_orbit");

	EarthMoonSun* EMS = new EarthMoonSun(jd);

	std::vector<double> jds;
	std::vector< std::array<double,3> > pos;
	
}

int main(int argc, char* argv[]) {

	int nSections = 30;
	double jd = 2460415.5;
	double Az = 10000;
	if(argc > 3) {
		nSections = (int)(std::stod(argv[1])/Section::SECTION_DAYS);
		jd = std::stod(argv[2]);
		Az = std::stod(argv[3]);
	}

	// const double jd0 = Util::getJDFromUTC(2024,4,15,0,0,0);
	test(nSections,jd,Az);
	
	return 0;
}
