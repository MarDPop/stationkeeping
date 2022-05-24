#include "../include/Simulation.h"
#include "../include/OrbitalDynamics.h"
#include "../include/OrbitalElements.h"
#include "../include/Util.h"
#include "../include/Config.h"
#include <cmath>
#include <array>

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

	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(cr3bp);
	ode.recording.set_record_interval(0.001);
	ode.set_timestep(1e-6);
	ode.run(x0,2);

	double period = tf*2/cr3bp->mean_motion;
	std::cout << "Period Guess: " << period/86400 << "days vs. " << cr3bp->get_period(Az)/86400 << "days." << std::endl;
	period = cr3bp->get_period(Az);

	Util::printOut(ode.recording,"output/CR3BP_orbit");

	Recording<6> orbit = OrbitComputation::get_cr3bp_halo_orbit(cr3bp,x0,600,period);

	Util::printOut(orbit,"output/CR3BP_orbit_2");
	
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
	test(nSections,jd,Az);
	
	return 0;
}
