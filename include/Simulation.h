#pragma once

#include <array>
#include <vector>
#include "OrbitalDynamics.h"
#include "Math.h"
#include "Matrix.h"
#include "ODE.h"
#include <thread>

struct Section {
	std::array<double,6> initial_state;
	double t_start, t_final;
	std::vector<double> times;
	std::vector< std::array<double,6> > states;
	Matrix<6,6> STM;
	ODE_RK4<6> ode;
	EarthMoonSun* dynamics;
	static const double SECTION_DAYS;
	
	Section(EarthMoonSun* f);
	
	void compute_states();

    double target(const std::array<double,6>& xp);

	void compute_STM();
};

class OrbitComputation {

	double julian_day;

public:

	OrbitComputation(){}
	~OrbitComputation(){}

	/**
	 * @brief Get the cr3bp halo object
	 * 
	 * @param dt 
	 * @param Az 
	 * @param sma 
	 * @param tol 
	 * @param mu1 
	 * @param mu2 
	 * @return std::array<double,6> 
	 */
	static void differential_correct_cr3bp(CR3BP* cr3bp, std::array<double,6>& x0, std::array<double,6>& xf, double& tf, bool top, const double& tol);

	static Recording<6> get_cr3bp_halo_orbit(CR3BP* cr3bp, std::array<double,6> x, double dt, double t_period);

	static void minimizeDX(std::vector<Section>& sections);

    static void minimizeDV(std::vector<Section>& sections);

	static void minimizeDV2(std::vector<Section>& sections);

	static void minimizeDV3(std::vector<Section>& sections);

    static void smooth(std::vector<Section>& sections);

    static double calcDV(std::vector<Section>& sections);

	/**
	 * @brief 
	 * 
	 * @param nSections 
	 */
	static void run_full_emphemeris(const int& nSections, const double& jd0);

	static void correct(std::vector<Section>& sections);


};

