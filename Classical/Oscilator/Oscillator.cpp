#pragma once

#include "Oscillator.h"

Oscillator::Oscillator(float mass, float k, float x0) {
	this->mass = mass;
	this->k = k;
	this->x0 = x0;
}

std::vector<pos_vel> Oscillator::RunOscillator(int tSteps, float dt, int mf, float distMax)
{

	std::vector<pos_vel> results(tSteps);


	if (mf == -1)
		mf = tSteps;


	float v = 0;
	float x = this->x0;

	for (int t = 0; t < tSteps; t++) {
		float acc = (-this->k * x) / this->mass; //: a = -kx/m
		v += acc * dt;

		if (t > 0 && t % mf == 0) {
			v += RandDisturbance(distMax);
		}
		
		
		x += v * dt;
		//Store position and velocity for this time step
		results[t] = pos_vel(x, v);
	}

	return results;
}

std::vector<pos_vel> Oscillator::RunMultipleOscillator(int tSteps, float dt, int mf, float distMax, int iterations, bool printProgress)
{
	std::vector<pos_vel> results(tSteps);
	float progress = 0;

	for (int i = 0; i < iterations; i++) {
		//Run a single iteration
		std::vector<pos_vel> it_res = this->RunOscillator(tSteps, dt, mf, distMax);
		//Add to results
		for (int t = 0; t < tSteps; t++) {
			results[t] += it_res[t];
		}


		if (printProgress) {
			progress = (i * 100.0) / iterations;
			std::cout << "(MF: " << mf << ", MaxD: " << distMax << ")Progress: " << int(progress) << "%\r";
			std::cout.flush();
		}

	}
	//Scale results by iterations count
	for (int t = 0; t < tSteps; t++) {
		results[t] /= (1.0 * iterations);
	}
	return results;
}

float Oscillator::CalculateAverageAmplitude(int tSteps, float dt, int iterations, int mf, float distMax, bool printProgress)
{

	float avgAmp = 0;


	float progress = 0;

	for (int it = 0; it < iterations; it++) {

		std::vector<pos_vel> results = this->RunOscillator(tSteps, dt, mf, distMax);
		avgAmp += CalculateAverageDisturbance(results);

		if (printProgress) {
			progress = (it * 100.0) / iterations;
			std::cout << "(MF: " << mf << ", MaxD: " << distMax << ")Progress: " << int(progress) << "%\r";
			std::cout.flush();
		}

	}
	avgAmp /= (1.0 * iterations);
	return avgAmp;
}

float RandDisturbance(float maxDist) {
	/// <summary>
	/// Returns a random value 'r':
	/// -maxDist <= r <= maxDist
	/// 
	/// </summary>
	/// <param name="maxDist">Maximum/minimum value to return</param>
	/// <returns>'r'</returns>
	return ((float(rand()) / RAND_MAX) - 0.5f) * 2 * maxDist;
}

float CalculateAverageDisturbance(std::vector<pos_vel> &pos_vel, float steadyStateFrac)
{

	int start = int(steadyStateFrac * pos_vel.size()) + 1;



	
	std::vector<float> turningPoints;

	for (int t = start; t < pos_vel.size(); t++) {
		//If the product of 2 differing time steps velocity is < 0, then we have a change in direction
		if (pos_vel[t - 1].vel * pos_vel[t].vel < 0) {
			turningPoints.push_back(0.5 * (pos_vel[t - 1].pos + pos_vel[t].pos));
		}
	}
	float sum = 0;
	for (int i = 0; i < turningPoints.size(); i++) {
		sum += abs(turningPoints[i]);
	}

	return sum / turningPoints.size();
}
