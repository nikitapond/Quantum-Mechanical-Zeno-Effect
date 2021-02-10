#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <thread>
#include "Oscillator.h"
#include <mutex>          // std::mutex


/// <summary>
/// Saves the vector of floats to the file path defined by 'save'
/// </summary>
/// <param name="res">Vector of floats to save</param>
/// <param name="save">File path to save to</param>
void SaveResults(std::vector<float> res, std::string save) {

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


const int T_STEPS = 640000;
const int ITERATIONS = 100;
const float dt = 0.001;


const int M_COUNT = 8;
const int DIST_COUNT = 21;
const int MeasureFreq[M_COUNT] = { -1, 64000, 32000, 16000, 8000, 4000, 2000, 1000 };
const float MaxDist[DIST_COUNT] = { 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
								   0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.40 };

//Array to hold all results before we save
float average_amplitudes[M_COUNT][DIST_COUNT];

int total_progress = 0;

/// <summary>
/// Starts a thread which will calculate the average amplitude for each different Max Disturbance (see MaxDist)
/// For each Max Disturbance, we run ITERATIONS number of times. Each iteration we measure the average amplitude
/// of the oscillator, We then take the average of this value for all iterations, and store in the array "a
/// </summary>
/// <param name="m"></param>
void StartAmpThread(int m) {
	int mf = MeasureFreq[m];


	for (int d = 0; d < DIST_COUNT; d++) {
		float di = MaxDist[d];
		Oscillator o(1, 1, 5);
		float avg = o.CalculateAverageAmplitude(T_STEPS, dt, ITERATIONS, mf, di, false);

		//Ensure thread safe access of array
		mtx.lock();
		total_progress++;
		int p_ = int(total_progress / (float(M_COUNT * DIST_COUNT)) * 100);

		std::cout << "Progress: " << p_ << "%\r";
		std::cout.flush();

		average_amplitudes[m][d] = avg;
		mtx.unlock();
	}


}

void RunAverageAmplitudeTest() {
	/// <summary>
	/// Starts the multithreading part of our simulations.
	/// Creates a new thread for each different frequency count (MeasureFreq), then runs for
	/// ITERATIONS number of times for each different disturbance (MaxDist)
	/// </summary>




	std::thread threads[M_COUNT];

	for (int m = 0; m < M_COUNT; m++) {
		threads[m] = std::thread(StartAmpThread, m);
	}
	//Wait for all threads to finish
	for (int i = 0; i < M_COUNT; i++) {
		threads[i].join();
	}
	std::cout << "results: " << std::endl;
	for (int m = 0; m < M_COUNT; m++) {
		for (int d = 0; d < DIST_COUNT; d++) {

			std::cout << average_amplitudes[m][d] << ",";

		}
		std::cout << "\n";
	}
	std::string save = "results/average_amplitudes2.csv";

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
				file << average_amplitudes[m][d] << "\n";
			else
				file << average_amplitudes[m][d] << ",";
		}
	}

	file.close();
	std::cout << "Succesfully saved results to " << save << std::endl;
}
