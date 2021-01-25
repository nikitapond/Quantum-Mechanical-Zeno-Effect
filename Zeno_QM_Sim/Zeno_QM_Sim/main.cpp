#pragma once
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <vector> 
#include <exception>
#include <mutex>
#include <fstream>
#include <math.h>       /* pow */
//#include <format>

#include <thread>
#include <chrono>

using namespace std;

class Timer {
private:
	std::chrono::steady_clock::time_point begin;
public:
	Timer() {
		begin = std::chrono::steady_clock::now();
	}
	~Timer() {
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

	}
};

float e = 2.7182818284;

string ToCVS(vector<float> res) {
	string str;

	for (int i = 0; i < res.size(); i++) {
		str.append(to_string(res[i]));
		if (i < res.size() - 1) {
			str.append(",");
		}
		
	}

	return str;
}
string BASE_PATH = "C:/Users/nikit/github/ZenoQM/Zeno_QM/Python/Results/";
void WriteString(string file, string s) {
	/*
	* Writes the string 's' into the file 'file'
	*/

	std::ofstream cvsFile;
	string path = BASE_PATH + file;
	
	cvsFile.open(path);
	cvsFile << s;
	cvsFile.close();
}

uint64_t shuffle_table[4] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
// The actual algorithm
inline uint64_t rand_fast(void)
{
	uint64_t s1 = shuffle_table[0];
	uint64_t s0 = shuffle_table[1];
	uint64_t result = s0 + s1;
	shuffle_table[0] = s0;
	s1 ^= s1 << 23;
	shuffle_table[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
	return result;
}

inline float very_randF() {
	return rand_fast() / double(UINT64_MAX);
}
inline float randF() {
	return rand() / double(RAND_MAX);
}

inline uint32_t bitExtracted(uint64_t number, int start, int len)
{
	unsigned mask;
	mask = ((1 << len) - start);
	return uint32_t(number & mask);
	//return (((1 << k) - 1) & (number >> (p - 1)));
}


vector<float> CalculateParticleDecay(int t_max, float dt, int particle_count, float Gamma) {
	/*
	We run a simulation, where we check each of our particles to see if, after a time passing of dt, 
	any given particle has decayed.
	The probability that a partical has decayed in a time period [t,t+dt] can be found by:
		p_decay(t, t+dt) = Gamma * dt * e^(-Gamma * dt)
	*/


	int length = int(t_max / dt);
	//In initial state, we have all particles

	int ten_perc = int(length * 0.1);
	vector<float> results(length + 1);
	results[0] = particle_count;
	bool needRan = true;
	uint64_t ran = 0;
	for (int i = 0; i < length; i++) {
		


		int decayed = 0;
		//Define the probability it will decay in this frame
		float decay_prob = dt * Gamma * pow(e, -Gamma * dt);

		for (int dec = 0; dec < results[i]; dec++) {

			
			if (rand_fast() / double(UINT64_MAX) < decay_prob) {
				decayed++;
			}
		}
		results[i+1] = results[i] - decayed;

		if (i % ten_perc == 0) {
			float cur_perc = int(100 *  (float(i) / length));
			cout << to_string(cur_perc) << "% complete, " << to_string(results[i + 1]) << " remaining" << endl;
		}

		if (results[i+1] <= 0)
			break;



	}
	return results;

}

void NormalParticleDecayTest() {
	/*
	This test shows that the rate of decay is independent of dt
	This 
	
	*/
	int num_particals = 1000000;
	float t_max = 1;

	{



		for (int t = 1; t < 5; t++) {
			//Scope based timer
			Timer time;
			float dt = 0.01 * pow(0.2, t);
			vector<float> res = CalculateParticleDecay(t_max, dt, num_particals, 10);
			res.insert(res.begin(), dt);
			string str = ToCVS(res);
			string file = "particle_decay";
			file.append(to_string(t));
			file.append(".csv");
			WriteString(file, str);
			cout << "Experiment " << to_string(dt) << "complete" << endl;
		}

	}
}

vector<float> HalfLifeDecay(float half_life, float t_max, float dt, int num_particles) {

	int length = int(t_max / dt);
	vector<float> results(length+1);

	results[0] = num_particles;
	for (int i = 1; i < length; i++) {
		int parts = int(num_particles * pow(0.5, (i * dt)/ half_life));
		results[i] = parts;
		if (parts < 0)
			break;

	}
	return results;

}

void HalfLifeDecayTest() {
	float half = 0.52084;
	int num_particals = 1000000;
	float t_max = 1;

	for (int t = 1; t < 5; t++) {
		//Scope based timer
		Timer time;
		float dt = 0.01 * pow(0.2, t);
		vector<float> res = HalfLifeDecay(half, t_max, dt, num_particals);
		res.insert(res.begin(), dt);
		string str = ToCVS(res);
		string file = "half_decay";
		file.append(to_string(t));
		file.append(".csv");
		WriteString(file, str);
		cout << "Experiment " << to_string(dt) << "complete" << endl;
	}

}
int main2() {

	NormalParticleDecayTest();
	HalfLifeDecayTest();
	
	return 0;
}