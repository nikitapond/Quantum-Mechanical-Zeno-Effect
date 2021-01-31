
#pragma once

#include <vector>
#include <string>
#include "RabiMCWF.h"
#include <fstream>
#include <iostream>
#include <iomanip>      // std::setprecision



void SaveResults(std::vector<Vec2R> res, std::string save) {



    std::ofstream file;
    file.open(save);

    for (int i = 0; i < res.size(); i++) {
        if (i == res.size() - 1) {
            file << res[i].x << "\n";
        }
        else {
            file << res[i].x << ",";
        }
    }
    for (int i = 0; i < res.size(); i++) {
        if (i == res.size() - 1) {
            file << res[i].y << "\n";
        }
        else {
            file << res[i].y << ",";
        }
    }
    file.close();
    std::cout << "Succesfully saved results to " << save << std::endl;
}

std::string roundStr(double val, int dp) {

    int mul = std::pow(10, dp);

    return std::to_string(roundf(val * mul) / mul);

}
std::string toStrMaxDecimals(double value, int decimals)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(decimals) << value;
    std::string s = ss.str();
    if (decimals > 0 && s[s.find_last_not_of('0')] == '.') {
        s.erase(s.size() - decimals + 1);
    }
    return s;
}
double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}

const double PI = 3.14159265359;
int main(){

    RabiMCWF system(1,10,1,2* PI);

    
    int total_iterations = 20000;
    int n_steps = 2000;
    int rabi_cycle_length = 500;

    double measurement_freq[] = { -1, 0.2, 0.1, 0.04, 16.0 / rabi_cycle_length,8.0 / rabi_cycle_length,4.0 / rabi_cycle_length, 2.0/ rabi_cycle_length };

    for (int i = 0; i < sizeof(measurement_freq) / sizeof(measurement_freq[0]); i++) {

        double mf = measurement_freq[i];
        int measureDelta;
        int iterations;

        std::string file_path = "mcwf_rabi_res_mf_";
        if (mf < 0){
            file_path += "inf.csv";
            measureDelta = n_steps; //Make no measurements
            iterations = 1;
        }        
        else {
            file_path += toStrMaxDecimals(mf,3) + ".csv";
            measureDelta = int(rabi_cycle_length * mf);
            iterations = total_iterations;
        }
            
       


        std::vector<Vec2R> res = system.RunExperiment(n_steps, measureDelta, 0, 0.001, iterations);
        std::cout << res.size() << std::endl;
        SaveResults(res, file_path);
    }

    return 0;
}
