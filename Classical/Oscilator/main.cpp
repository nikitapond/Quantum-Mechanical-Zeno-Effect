//#include <vector>
#include <fstream>
#include <string>
#include <iostream>
//#include <thread>
//#include "Oscillator.h"
//#include <mutex>          // std::mutex
#include <filesystem>



//Function in AmplitudeCalculation.cpp that starts a multi threaded 
//calculation
void RunAverageAmplitudeTest();
void RunPositionVelocityCalc();


const int T_STEPS = 640000;
const int ITERATIONS = 100;
const float dt = 0.001;


const int M_COUNT = 8;
const int DIST_COUNT = 21;
const int MeasureFreq[M_COUNT] = { -1, 6400, 3200, 1600, 800, 400, 200, 100 };
const float MaxDist[DIST_COUNT] = { 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
								   0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.40 };


void writeAllConfigurations(std::string save) {
	std::ofstream file;
	file.open(save);


	file.close();


}


int main() {

	namespace fs = std::filesystem;
	fs::create_directory("results");
	fs::create_directory("results/pos_vel");


	//RunAverageAmplitudeTest();
	RunPositionVelocityCalc();
	
	


	return 0;
};