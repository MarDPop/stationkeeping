#pragma once

#include <array>
#include <vector>
#include "OrbitalDynamics.h"
#include "Math.h"
#include "Matrix.h"
#include "Ode.h"
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

    void target(const std::array<double,6>& xp);

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
	 * @return std::vector< std::array<double,6> > 
	 */
	static const std::vector< std::array<double,6> > get_cr3bp_halo(const double& dt, const double& Az, const double& sma, const double& tol, const double& mu1, const double& mu2);

	static void minimizeDX(std::vector<Section>& sections);

    static void minimizeDV(std::vector<Section>& sections);

    static void smooth(std::vector<Section>& sections);

    static double calcDV(std::vector<Section>& sections);

	/**
	 * @brief 
	 * 
	 * @param nSections 
	 */
	static void run_full_emphemeris(const int& nSections);

	/**
	 * @brief 
	 * 
	 * @param cr3bp 
	 * @param dynamics 
	 * @param state_guess 
	 * @param jd 
	 * @return std::array<double,6> 
	 */
	static std::array<double,6> convert(CR3BP* cr3bp, EarthMoonSun* dynamics, const std::array<double,6>& state_guess, const double& jd);

};

