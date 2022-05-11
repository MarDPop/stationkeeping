#pragma once

#include <string>
#include <vector>
#include <array>

class OrbitalElements {
public:
	const double mu;
	std::vector< double > jds;
	std::vector< std::array<double,6> > oes;
	std::vector< std::array<double,6> > d_oe;
	std::vector< double > mean_anomaly;
	std::vector< double > mean_motion;
	int idx;

	static const int LIGHTSPEED; 
	static const double EARTH_MU; 
	static const double MOON_MU; 
	static const double SUN_MU; 
	static const double JUPITER_SYSTEM_MU;
	static const double MARS_SYSTEM_MU;
	static const double ICRF_MU;
	static const double EMBARY_MU;
	static const double SOLAR_IRRADIANCE_AU; 
	static const double SOLAR_PRESSURE_AU;
	static const double AU; 
	static const double UT_TDB_DELTAT;
	
	OrbitalElements( const std::vector<double>& JD, const std::vector< std::array<double,6> >& OE, const double& MU);
	OrbitalElements(const std::string& fn, const double& MU);
	~OrbitalElements();
	
	inline double getMu() const{
		return this->mu;
	}
	
	static double trueAnomalyFromEccentricAnomaly( const double& EA, const double& eccentricity );

    static double trueAnomalyFromMeanAnomaly( const double& MA, const double& eccentricity );
	
	static double trueAnomalyFromMeanAnomalyApprox( const double& MA, const double& eccentricity );
    
    static double meanAnomalyFromTrueAnomalyApprox( const double& f, const double& eccentricity);

    static double eccentricAnomalyFromMeanAnomaly( const double& MA, const double& eccentricity);
    
    static double meanAnomalyFromEccentricAnomaly( const double& EA, const double& eccentricity);
    
    static double eccentricAnomalyFromTrueAnomaly( const double& TA, const double& eccentricity );
    
    static double meanAnomalyFromTrueAnomaly( const double& f, const double& eccentricity);
    
    static std::array<double,6> kepler2cartesian(const std::array<double,6>& oe, const double& mu);
	
	static std::array<double,3> kepler2position(const std::array<double,6>& oe);
    
    static std::array<double,6> cartesian2kepler(const std::array<double,6>& state, const double& mu);
	
	std::array<double,6> get(const double& jd) const;
	
	std::array<double,6> getSimple(const double& jd);
};