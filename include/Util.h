#pragma once

#include <sstream>

namespace Util {
	
	constexpr int J2000 = 2451545;
	constexpr int JULIAN_DAY = 86400;
	constexpr double SIDREAL_DAY = 86164.09054;
	
	void writeCSV(const std::vector< std::vector<double> > & data, const std::string & filename) {
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
    
    std::vector<std::string> split(const std::string& s, const char& delimiter) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }
	
	int getJDNFromGregorianYMD(int y, int m, int d){
        int m12 = (m - 14)/12;
        return  (1461*(y + 4800 + m12))/4 +
                (367*( m - 2 - 12 * m12 ))/12 -
                (3*( (y + 4900 + m12)/100))/4 +
                d - 32075;
    }
	
	double getJDFromUTC(int y, int m, int d, int hour, int min, double sec){
        int jdn = getJDNFromGregorianYMD(y,m,d);
		double day_sec = (hour - 12)*3600 + min*60 + sec;
		return jdn + day_sec/86400.0;
    }
	
}