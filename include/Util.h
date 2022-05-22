#pragma once

#include "Dynamics.h"
#include "ODE.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

namespace Util {
	
	constexpr int J2000 = 2451545;
	constexpr int JULIAN_DAY = 86400;
	constexpr double SIDREAL_DAY = 86164.09054;
	
	inline void writeCSV(const std::vector< std::vector<double> > & data, const std::string & filename) {
        std::ofstream myfile;
        myfile.open (filename);
        for(std::size_t i = 0; i < data.size();i++){
            std::vector<double> row = data[i];
            myfile << row[0];
            for(std::size_t j = 1; j < row.size();j++){
                myfile << "," << row[j];
            }
            myfile << "\n";
        }
        myfile.close();
    }
    
    inline std::vector<std::string> split(const std::string& s, const char& delimiter) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }

    inline void printOut(const Recording<6>&  record, std::string filename){
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

    inline void printOut(const std::vector<double>& t,const std::vector<std::array<double,3> >& x, std::string filename){
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

    inline void printOut(EarthMoonSun* dynamics, const std::vector<Section>& sections, std::string filename){
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
	
	inline int getJDNFromGregorianYMD(int y, int m, int d){
        int m12 = (m - 14)/12;
        return  (1461*(y + 4800 + m12))/4 +
                (367*( m - 2 - 12 * m12 ))/12 -
                (3*( (y + 4900 + m12)/100))/4 +
                d - 32075;
    }
	
	inline double getJDFromUTC(int y, int m, int d, int hour, int min, double sec){
        int jdn = getJDNFromGregorianYMD(y,m,d);
		double day_sec = (hour - 12)*3600 + min*60 + sec;
		return jdn + day_sec/86400.0;
    }
	
}