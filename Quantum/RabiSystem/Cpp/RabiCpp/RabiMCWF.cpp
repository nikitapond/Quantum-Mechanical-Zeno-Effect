#pragma once

#include "RabiMCWF.h"
#include "Vec2.h"
#include "Mat2.h"
#include <random>

// Returns a number between 0 (inclusive) and 1 (exclusive)
double MyRand(){
    return ((double)rand())/(RAND_MAX+1);
}
/// <summary>
/// Define various matricies and operators used in our calculations
/// </summary>
Mat2 Sp(0, 0, 1, 0); //Sigma+
Mat2 Sm(0, 1, 0, 0); //Sigma-
Mat2 Ee(0, 0, 0, 1); //S+ * S-
Mat2 Ve(0, 1, 1, 0); //S+ + S-
Mat2 Id(1, 0, 0, 1); //Identity matrix

Mat2 aDag(0, 0, 1, 0);
Mat2 a(0, 1, 0, 0);

Mat2 aDagMinusA(0, -1, 1, 0); 
Mat2 aDagMultA(0, 0, 0, 1);
/// <summary>
/// Default constructor
/// </summary>
RabiMCWF::RabiMCWF() {
    hbar = 1.;
    n_th = 1.;
    kappa = 1.;
    eta = 1.;
    Init();
}
/// <summary>
/// Main constructor
/// </summary>
/// <param name="hbar">value of hbar to use. Best set to 1</param>
/// <param name="n_th">Abstract value used in defining magnitude of jump operators</param>
/// <param name="kappa">Constant relating to time, related to magnitude of jump operators.
/// Best set to 1.</param>
/// <param name="eta">Rabi frequency of system</param>
RabiMCWF::RabiMCWF(double hbar, double n_th, double kappa, double eta){
    this->hbar = hbar;
    this->n_th = n_th;
    this->kappa = kappa;
    this->eta = eta;
    Init();
}
/// Calls all functions required to initialise our system for calculations.
void RabiMCWF::Init() {
    this->ham = GenerateHamiltonian();
    this->jumpOperators = GenerateJumpOperators();
    this->hamNH = GenerateNonHermitianHamiltonian();
}

//Generates the hamiltonian for this system
Mat2 RabiMCWF::GenerateHamiltonian(){
    return aDagMinusA * (this->eta * I);
}

//Generates the non-hermitian hamiltonian for this system
Mat2 RabiMCWF::GenerateNonHermitianHamiltonian(){
    Mat2 jump_part(0,0,0,0);

    for(int i=0; i<2; i++){
        Mat2 jump_H = (this->jumpOperators[i]).Hermitian();
        jump_part += DOT(jump_H, this->jumpOperators[i]);
    }
    return this->ham - jump_part * ( 0.5 * this->hbar * I);
}
// Generates the jump operators for this system
std::vector<Mat2> RabiMCWF::GenerateJumpOperators(){
    std::vector<Mat2> jumps(2);

    Mat2 j_up = a * sqrt(2 * this->kappa * (this->n_th));
    Mat2 j_down = aDag * sqrt(2 * this->kappa * (this->n_th + 1));

    jumps[0] = j_up;
    jumps[1] = j_down;
    return jumps;
}
/// <summary>
/// Evolves the system according to a quantum jump caused by an applied measurement
/// </summary>
/// <param name="dt">Time step of simulatuon</param>
/// <param name="dp">The relative probability of this jump occuring</param>
/// <param name="it">Current iteration</param>
/// <param name="phi">Reference to the vector of time evolution of phi </param>
void RabiMCWF::NoJump(double dt, double dp, int it, std::vector<Vec2>& phi){
    Mat2 lhs = Id - (this->hamNH * (I * dt/this->hbar));

    Vec2 res = lhs.dot(phi[it-1]);
    res /= sqrt(1-dp);
    phi[it] = res;
}
/// <summary>
/// Evolves our system according to the Hamiltonian, in the case that no jump occurs
/// </summary>
/// <param name="dt"></param>
/// <param name="p_k"></param>
/// <param name="jumpOp"></param>
/// <param name="it"></param>
/// <param name="phi"></param>
void RabiMCWF::Jump(double dt, double p_k, Mat2& jumpOp, int it, std::vector<Vec2>& phi){
    double mult = sqrt(dt/p_k);

    Vec2 rhs = (jumpOp.dot(phi[it-1]));

    phi[it] = rhs * mult;
}
/// <summary>
/// Calculates the probability of each of the two jump operators for the current
/// quantum superposition.
/// </summary>
/// <param name="phi">Current quantum superposition</param>
/// <returns>Vec2R.x -> probability of jump[1],
///          Vec2R.y -> probability of jump[0]</returns>
Vec2R RabiMCWF::CalcPk(Vec2& phi){    
    double y = real(phi.Conj().dot(this->jumpOperators[0].Hermitian()).dot(this->jumpOperators[0]).dot(phi));
    double x = real(phi.Conj().dot(this->jumpOperators[1].Hermitian()).dot(this->jumpOperators[1]).dot(phi));

    return Vec2R(x,y);
}

/// <summary>
/// Main loop for a single given evolution. We start our quantum system in a known state (ground state) 
/// and allow it to evolve over time, taking measurements (collapsing quantum state) when
/// required.
/// </summary>
/// <param name="dt">Time step for each iteration</param>
/// <param name="nSteps">Total number of time steps to simulate</param>
/// <param name="phi">Reference to a vector of Vec2 that will store results</param>
/// <param name="measureDelta">Time steps between measurements</param>
/// <returns>True if sucessfull, false if error occurs</returns>
bool RabiMCWF::Propagate(double dt, int nSteps, std::vector<Vec2>& phi, int measureDelta){
    phi[0].x = 1;
    phi[0].y = 0;
    int last_delta = 0;

    for(int it=1; it<nSteps; it++){
   
        //Calculate the relative probability of each jump
        Vec2R pK_all = this->CalcPk(phi[it-1]) * dt;

        double dp = pK_all.x + pK_all.y;
        //Generate random number for 
        double epsilon = MyRand();
        
        bool measure=false;

        if(it - last_delta >= measureDelta){
            measure = true;
            last_delta = it;

        }
        
        if(measure){
            double p_current = 0;
            pK_all /= dp;
            double r = MyRand();

            bool fin= false;
            for(int i=0; i<2; i++){
                if(fin){ break; }

                double pk = pK_all[i];
                p_current += pk;
               // std::cout << "all rand: " << pK_all.ToString() << "\n," << pk << "," << p_current << ", " << r << std::endl;
                if(p_current >= r){
                    this->Jump(dt, pk, this->jumpOperators[i], it, phi);
                    phi[it].normalize();
                    //std::cout << "After a jump, our result was " << phi[it].ToString() << std::endl;
                   // std::cout << "norm: " << phi[it].ToString() << std::endl;
                    fin=true;
                    break;
                }
            }
            //Used in bebuging, this should never happen
            if (!fin) {

                std::cout << "Jump has failed!\n  r: " << r << "\nPk_all: " << pK_all.ToString()  << "\n r: " << r << std::endl;
                this->NoJump(dt, dp, it, phi);

                phi[it].normalize();
                return false;
            }
            //only apply measurements for a single time step
            measure = false;

        }else{
            //If no measurement happens, evolve the system normally
            this->NoJump(dt, dp, it, phi);
            phi[it].normalize();
   
        }

    }

}


std::vector<Vec2R> RabiMCWF::RunExperiment(int nSteps, int measureDelta, double dt, int iterations){

    //Vector to store the final results (probability of being in excited state at given time step)
    std::vector<Vec2R> results(nSteps);
    //Vector to store the wave function result for each iteration
    std::vector<Vec2> phi(nSteps);
    //Set random seed for repeatability
    srand(time(NULL));
    //Run all iterations
    float progress = 0;
    for(int is=0; is<iterations; is++){

        //calculate all our phi
        if (!this->Propagate(dt, nSteps, phi, measureDelta))
            return results;
        //Store the probability (|psi|^2) at each iteration
        for(int it=0; it<nSteps; it++){
            results[it].x += (real(conj(phi[it].x) * phi[it].x));
            results[it].y += (real(conj(phi[it].y) * phi[it].y));
        }
        progress = (is * 1.0) / iterations;
        std::cout << "Progress: " << int(progress * 100) << "%\r";
        std::cout.flush();
    }
    //Scale results
    for (int it = 0; it < nSteps; it++) {
        results[it].x /= (1.0 * iterations);
        results[it].y /= (1.0 * iterations);
    }
    
    return results;
}