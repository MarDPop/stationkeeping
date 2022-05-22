#include "../include/Simulation.h"


void test(const std::size_t& nSections) {
	Simulation sim;
	sim.run(nSections);
	
}

int main(int argc, char* argv[]) {

	int nSections = 30;
	if(argc > 1) {
		nSections = (int)(std::stod(argv[1])/Section::SECTION_DAYS);
	}
	test(nSections);
	
	return 0;
}
