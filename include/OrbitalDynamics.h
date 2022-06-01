#pragma once

#include "Dynamics.h"

#include "Matrix.h"
#include <array>
#include <vector>
#include <string>

class EarthMoonSun;

class CR3BP;

class OrbitalDynamics : public Dynamics<6> {

public: 
	virtual void getA(const std::array<double,6>& x, const double& jd, double** A) const {}

	static std::array<double,3> convert_cr3bp_to_inertial_pos(EarthMoonSun* dynamics, const std::array<double,6>& state, const double& jd);

	static std::array<double,3> convert_cr3bp_to_rotating_barycenter(EarthMoonSun* dynamics, const std::array<double,6>& state, const double& jd);

	/**
	 * @brief 
	 * 
	 * @param cr3bp 
	 * @param dynamics 
	 * @param state_guess 
	 * @param jd 
	 * @return std::array<double,6> 
	 */
	static std::array<double,6> convert_cr3bp_to_inertial(EarthMoonSun* dynamics, const std::array<double,6>& state_guess, const double& jd);

};

class CR3BP : public OrbitalDynamics {
	
public:

	const double mu;
	const double mu1;
	const double mu_body1;
	const double mu_body2;
	const double sma;
	const double mean_motion;

	static constexpr double EARTH_L1 = 0.836915155115;
	
	CR3BP();
	CR3BP(const double& MU1, const double& MU2, const double& SMA);
	
	void getA(const std::array<double,6>& x, const double& jd, double** A) const;
	
	double getCL1(const double& yL,const int& n);

	double getCL2(const double& yL,const int& n);
	
	double getCL3(const double& yL,const int& n);

	double getL1();
	
	double getL2();
	
	double getL3();

	double get_period(const double& Az0);
	
	std::array<double,6> get_halo_initial_state_3rd_order(const double& Az0, const double& phi, const double& time, const int& n);
	
	std::array<double,6> convert_state_to_inertial(const std::array<double,6>& x) ;

	std::array<double,6> get_state_rate(const std::array<double,6>& x, const double& t) const;

};

struct VectorTable {
	std::vector< double > jds;
	std::vector< std::array< double,3 > > pos;
	std::vector< std::array< double,3 > > vel;
	
	VectorTable(const std::string& fn);
	
	std::array<double,3> getPos(const double& jd) const;

	std::array<double,3> getVel(const double& jd) const;
	
};

class EarthMoonSun : public OrbitalDynamics {
	
	std::vector<double> jds;
	std::vector<std::array<double,3> > r_earth;
	std::vector<std::array<double,3> > r_moon;
	std::vector<std::array<double,3> > r_sun;
	std::vector<std::array<double,3> > r_emb;
	std::vector<std::array<double,3> > v_emb;
	std::vector<std::array<double,3> > omega;
	
	std::vector<std::array<double,3> > d_earth;
	std::vector<std::array<double,3> > d_moon;
	std::vector<std::array<double,3> > d_sun;
	std::vector<std::array<double,3> > d_emb;
	std::vector<std::array<double,3> > dv_emb;
	
	const double djd = 0.1;
	const double jd0;
	double jdf;
	
public:	
	const VectorTable* earth;
	const VectorTable* moon;
	const VectorTable* emb;
	const VectorTable* sun;
	
	EarthMoonSun(const double& JD0);
	
	inline double getJD0() const{
		return this->jd0;
	}

	inline double getJDf() const{
		return this->jdf;
	}
	
	std::array< std::array<double,3>, 4> getEarthMoonBarycenterCS(const double& jd) const;

	std::array< std::array<double,3>, 4> getEarthMoonL1CS(const double& jd) const;
	
	std::array<double,6> cr3bp_to_embj2000(const std::array<double,6>& cr3bp, const double& jd) const;
	
	void getA(const std::array<double,6>& x, const double& jd, double** A) const;
	
	void getA(const std::array<double,6>& x, const double& jd, Matrix<6,6>& A) const;
	
	Matrix<3,3> getG(const std::array<double,6>& x, const double& jd) const ;

	std::array<double,6> get_state_rate(const std::array<double,6>& x, const double& t) const;
};
