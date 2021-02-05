#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <thread>
#include "Oscillator.h"
#include <mutex>          // std::mutex

void SaveResults(std::vector<float> res, std::string save){

	std::ofstream file;
	file.open(save);
	for (int i = 0; i < res.size(); i++) {
		if (i == res.size() - 1)
			file << res[i] << "\n";
		else
			file << res[i] << ",";
	}
	file.close();
	std::cout << "Succesfully saved results to " << save << std::endl;
}


std::mutex mtx;           // mutex for critical section


const int T_STEPS = 64000;
const int ITERATIONS = 500;
const float dt = 0.001;

const int M_COUNT = 8;
const int DIST_COUNT = 21;
const int MeasureFreq[M_COUNT] = { -1, 6400, 3200, 1600, 800, 400, 200, 100 };
const float MaxDist[DIST_COUNT] = { 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
								   0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.40 };

float all_results[M_COUNT][DIST_COUNT];


void StartThread(int m) {
	int mf = MeasureFreq[m];


	for (int d = 0; d < DIST_COUNT; d++) {
		float di = MaxDist[d];
		Oscillator o(1, 1, 5);
		float avg = o.CalculateAverageAmplitude(T_STEPS, dt, ITERATIONS, mf, di, false);
		mtx.lock();

		all_results[m][d] = avg;
		std::cout << "Thread: " << mf << "," << di << " complete: " << avg << std::endl;
		mtx.unlock();
	}

	
}

int main() {



	
	int t_count = 8;

	std::thread threads[M_COUNT];

	int t = 0;
	for (int m = 0; m < M_COUNT; m++) {
		threads[t] = std::thread(StartThread, m);
		t++;
	}

	for (int i = 0; i < t; i++) {
		threads[i].join();
	}
	std::cout << "results: " << std::endl;
	for (int m = 0; m < M_COUNT; m++) {
		for (int d = 0; d < DIST_COUNT; d++) {

			std::cout << all_results[m][d] << ",";
			
		}
		std::cout << "\n";
	}
	std::string save = "all_results.csv";

	std::ofstream file;
	file.open(save);
	//We print out in a way that pd data frame can easily read.
	//The top left (first) cell is empty

	file << "MFreqs,";
	//We do the headers like this
	for (int d = 0; d < DIST_COUNT; d++) {
		if (d == DIST_COUNT - 1)
			file << MaxDist[d] << "\n";
		else
			file << MaxDist[d] << ",";
	}
	for (int m = 0; m < M_COUNT; m++) {
		file << MeasureFreq[m] << ",";
		for (int d = 0; d < DIST_COUNT; d++) {

			if (d == DIST_COUNT - 1)
				file << all_results[m][d] << "\n";
			else
				file << all_results[m][d] << ",";
		}
	}

	file.close();
	std::cout << "Succesfully saved results to " << save << std::endl;


	/*
	int measureFrq[] = { -1, 6400, 3200, 1600, 800, 400, 200, 100 };

	std::vector<float> results(8);

	Oscillator o(1, 1, 5);


	std::vector<std::thread> threads;

	int tSteps = 640000; //640 thousand
	int iterations = 1000;
	for (int i = 0; i < 8; i++) {
		int mf = measureFrq[i];
		float avg = o.CalculateAverageAmplitude(tSteps, 0.001, iterations, mf, 0.25);
		results[i] = avg;
		std::cout << "average amplitude was " << avg << " for measurement frequency: " << mf << std::endl;

	}
	SaveResults(results, "AvgAmp_const_dist.csv");
	*/
	return 0;
};