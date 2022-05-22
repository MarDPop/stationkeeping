#include "../include/Simulation.h"

void test(const int& nSections) {
	OrbitComputation::run_full_emphemeris(nSections);
	
}

int main(int argc, char* argv[]) {

	int nSections = 30;
	if(argc > 1) {
		nSections = (int)(std::stod(argv[1])/Section::SECTION_DAYS);
	}
	test(nSections);
	
	return 0;
}
