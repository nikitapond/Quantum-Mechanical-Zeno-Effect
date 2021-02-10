#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <thread>
#include "Oscillator.h"
#include <mutex>          // std::mutex
#include <cstdlib>
//#include <experimental/filesystem>
#include <filesystem>
std::mutex mtx_1;           // mutex for critical section


//Array to hold all results before we save
const int T_STEPS = 640000;
const int ITERATIONS = 2500;
const float dt = 0.001;


const int M_COUNT = 8;
const int DIST_COUNT = 21;
const int MeasureFreq[M_COUNT] = { -1, 64000, 32000, 16000, 8000, 4000, 2000, 1000 };
const float MaxDist[DIST_COUNT] = { 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
								   0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.40 };



void SaveResult(std::vector<pos_vel> res, std::string save) {
	std::ofstream file;
	file.open(save);

	file << "pos,vel\n"; //Collumn titles

	for (int t = 0; t < res.size(); t++) {
		file << res[t].pos << "," << res[t].vel << "\n";
	}
	/*
	for (int t = 0; t < res.size(); t++) {
		if (t == res.size() - 1)
			file << res[t].pos << "\n";
		else
			file << res[t].pos << ",";
	}

	file << "vel, ";
	for (int t = 0; t < res.size(); t++) {
		if (t == res.size() - 1)
			file << res[t].vel << "\n";
		else
			file << res[t].vel << ",";
	}*/

	file.close();
	std::cout << "Succesfully saved results to " << save << std::endl;

}


float p = 0;

void StartPVThread(int m) {


	int d = 20;
	Oscillator o(64, 1, 5); 
	//omega = sqrt(k/m) = 0.001
	//dt = 0.001 -> 1 oscillation every 
	std::vector<pos_vel> pos_vel = o.RunMultipleOscillator(T_STEPS, dt, MeasureFreq[m], MaxDist[d], 1, false);
	std::string fName = "results/pos_vel/mf_" + std::to_string(MeasureFreq[m]) + "_md_" + std::to_string(MaxDist[d]) + "_.csv";
	SaveResult(pos_vel, fName);
	/*
	for (int d = 0; d < DIST_COUNT; d++) {
		Oscillator o(10, 0.1, 5);
		std::vector<pos_vel> pos_vel = o.RunMultipleOscillator(T_STEPS, dt, MeasureFreq[m], MaxDist[d], 1, false);
		std::string fName = "results/pos_vel/mf_" + std::to_string(MeasureFreq[m]) + "_md_" + std::to_string(MaxDist[d]) + "_.csv";
		SaveResult(pos_vel, fName);
		mtx_1.lock();

		p++;
		std::cout << "Progress: " << int(100 *p / (M_COUNT * DIST_COUNT)) << "%\r";
		std::cout.flush();
		mtx_1.unlock();

	}*/
}

void RunPositionVelocityCalc() {

	std::thread threads[M_COUNT];

	for (int m = 0; m < M_COUNT; m++) {
		threads[m] = std::thread(StartPVThread, m);
	}
	//Wait for all threads to finish
	for (int i = 0; i < M_COUNT; i++) {
		threads[i].join();
	}
	std::cout << "complete" << std::endl;




	

}