
#pragma once

#include <vector>
#include <string>
#include "RabiMCWF.h"
#include <fstream>
#include <iostream>
#include <iomanip>      // std::setprecision

const double PI = 3.14159265359;

/// <summary>
/// Saves the results of our simulation to the specified save file
/// as a csv file
/// </summary>
/// <param name="res">Vector of real vec2 data-points</param>
/// <param name="save">The file path to save the results to</param>
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
/// <summary>
/// Rounds the value to the specified number of decimal places,
/// then returns this value as a string.
/// </summary>
/// <param name="value">Value we wish to round</param>
/// <param name="decimals">Number of dp we wish</param>
/// <returns>'value' rounded to decimals' dp, as a string</returns>
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

/// <summary>
/// Starts our calculations
/// </summary>
/// <returns></returns>
int main(){

    //Define system
    RabiMCWF system(1,10,1,2* PI);

    
    int total_iterations = 20000; //Number of times to repeat and take an average
    int n_steps = 2000; //Number of time steps to evolve the system over
    int rabi_cycle_length = 500; //2000 timesteps = 5 rabi cycles, 500 = 1 rabi cycle
    //Define number of iterations between each measurement -> smaller value = mroe measurements
    double measurement_freq[] = { -1, 0.2, 0.1, 0.04, 16.0 / rabi_cycle_length,8.0 / rabi_cycle_length,4.0 / rabi_cycle_length, 2.0/ rabi_cycle_length };

    for (int i = 0; i < sizeof(measurement_freq) / sizeof(measurement_freq[0]); i++) {

        double mf = measurement_freq[i];
        int measureDelta;
        int iterations;
        //Calculated needed values
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
        //Run experiment and store results
        std::vector<Vec2R> res = system.RunExperiment(n_steps, measureDelta, 0.001, iterations);
        std::cout << res.size() << std::endl;
        SaveResults(res, file_path);
    }

    return 0;
}
