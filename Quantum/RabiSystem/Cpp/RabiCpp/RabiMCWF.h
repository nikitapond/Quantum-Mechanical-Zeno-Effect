#pragma once

#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "Mat2.h"

const int N = 2;




struct Vec2R{
    double x, y;
    Vec2R(){
        x=0;
        y=0;
    }
    Vec2R(double x, double y){
        this->x = x;
        this->y = y;
    }
    Vec2R operator*(double t){
        return Vec2R(x*t, y*t);
    }
    Vec2R operator/=(double t){
        this-> x/= t;
        this -> y /=t;
        return *this;
    }
    double operator[](int index){
        if(index == 0) return this->x;
        if(index == 1) return this->y;
        return 0; 
    }
    std::string ToString() {
        return std::to_string(x) + ", " + std::to_string(y);
    }
};

class RabiMCWF{

private:
    double hbar, n_th, kappa, eta;

    Mat2 ham; //Hamiltonian for this system
    Mat2 hamNH; //Non hermitian hamiltonian
    std::vector<Mat2> jumpOperators; //Jump operators

    void Init();

    Mat2 GenerateHamiltonian();
    Mat2 GenerateNonHermitianHamiltonian();
    std::vector<Mat2> GenerateJumpOperators();

    void NoJump(double dt, double dp, int it, std::vector<Vec2>& phi);
    void Jump(double dt, double p_k, Mat2& jumpOp, int it, std::vector<Vec2>& phi);

    Vec2R CalcPk(Vec2& phi);

    bool Propagate(double dt, int nSteps, std::vector<Vec2>& phi, int measureDelta);

public:
    RabiMCWF();
    RabiMCWF(double hbar, double n_th, double kappa, double eta);

    std::vector<Vec2R> RunExperiment(int nSteps, int measureDelta, double dt, int iterations);

};
