#include "../include/Simulation.h"
#include "../include/OrbitalDynamics.h"
#include "../include/OrbitalElements.h"
#include "../include/Util.h"
#include "../include/Config.h"
#include "../include/Math.h"
#include <cmath>
#include <array>
#include <iomanip>

void printSections(const std::vector<Section>& sections, std::string suffix){

	if(sections.size() < 1){
		return;
	}
	
	EarthMoonSun* dynamics = sections[0].dynamics;
	const double jd0 = dynamics->getJD0();
	for(int i = 0; i < sections.size(); i++){
		const Recording<6>& recording = sections[i].ode.recording;
		Util::printOut(recording,"output/section_inertial_" + std::to_string(i) + "_" + suffix);
		std::vector< double > jds;
		std::vector< std::array<double,3> > rotating_positions;

		for(int j = 0; j < recording.number_entries();j++){
			double t = recording.time_at(j);
			std::array<double,6> x = recording.state_at(j);
			double jd = t/Util::JULIAN_DAY + jd0;
			std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonL1CS(jd);
			std::array<double,3> r = {x[0] - frame[0][0], x[1] - frame[0][1], x[2] - frame[0][2]};
			std::array<double,3> r_L1;
			r_L1[0] = Math::dot(frame[1],r);
			r_L1[1] = Math::dot(frame[2],r);
			r_L1[2] = Math::dot(frame[3],r);

			jds.push_back(jd);
			rotating_positions.push_back(r_L1);
		}

		Util::printOut(jds,rotating_positions,"output/section_rotating_" + std::to_string(i) + "_" +suffix);
	}
}


void test(const int& nSections, const double& jd, const double& Az) {
	VectorTable earth = VectorTable("resources/EARTH_EMB_VECTOR.dat");
    VectorTable moon = VectorTable("resources/MOON_EMB_VECTOR.dat");

	double sma = 381000; // approximately average

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

	std::cout << "Printing converted CR3BP Orbits." << std::endl;

	std::vector< double > jds;
	std::vector< std::array<double,3> > positions;
	std::vector< std::array<double,3> > BC_positions;
	double time = 0;
	double dt = 60;
	const double time_final = nSections*Section::SECTION_DAYS*Util::JULIAN_DAY;
	while(time < time_final){
		double t = fmod(time, period) * cr3bp->mean_motion;
		double jdi = jd + time/Util::JULIAN_DAY;
		positions.push_back(OrbitalDynamics::convert_cr3bp_to_inertial_pos(EMS,orbit.get(t),jdi));
		BC_positions.push_back(OrbitalDynamics::convert_cr3bp_to_rotating_barycenter(EMS,orbit.get(t),jdi));
		jds.push_back(jdi);
		time += dt;
	}

	std::vector< std::array<double,3> > earth_pos;
	std::vector< std::array<double,3> > moon_pos;
	for(const double& jdi : jds){
		earth_pos.push_back(EMS->earth->getPos(jdi));
		moon_pos.push_back(EMS->moon->getPos(jdi));
	}
	
	Util::printOut(jds,positions,"output/CR3BP_orbit_converted");
	Util::printOut(jds,earth_pos,"output/earth_orbit");
	Util::printOut(jds,moon_pos,"output/moon_orbit");
	Util::printOut(jds,BC_positions,"output/CR3BP_orbit_converted_BC");

	std::vector< Section > sections(nSections,EMS);

	time = 0;
	dt = Section::SECTION_DAYS*86400;
	for(int i = 0; i < nSections; i++){
		double t = fmod(time, period) * cr3bp->mean_motion;
		double jdi = jd + time/Util::JULIAN_DAY;
		sections[i].initial_state = OrbitalDynamics::convert_cr3bp_to_inertial(EMS,orbit.get(t),jdi);
	
		sections[i].t_start = time;
		time += dt;
		sections[i].t_final = time;
		sections[i].compute_states();
	}

	std::cout << "Printing initial sections." << std::endl;
	printSections(sections,"orig");

	OrbitComputation::minimizeDX(sections);

	printSections(sections,"minDx");

	try {

		OrbitComputation::minimizeDV3(sections);

	} catch ( const char* e) {
		std::cerr << e << std::endl;
		return;
	}

	printSections(sections,"minDv");

	OrbitComputation::minimizeDX(sections);

	for ( int i = 0; i < 5; i++) {
		OrbitComputation::minimizeDV3(sections);

		OrbitComputation::minimizeDX(sections);
	}

	printSections(sections,"pass_final");

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
