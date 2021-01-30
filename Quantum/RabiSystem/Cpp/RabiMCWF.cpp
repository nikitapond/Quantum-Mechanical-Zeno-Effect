#pragma once

#include "RabiMCWF.h"

double MyRand(){
    return std::rand()/RAND_MAX;
}


RabiMCWF::RabiMCWF(double hbar, double n_th, double kappa, double eta){
    this->hbar = hbar;
    this->n_th = n_th;
    this->kappa = kappa;
    this->eta = eta;
    this->ham = GenerateHamiltonian();
    this->jumpOperators = GenerateJumpOperators();
    this->hamNH = GenerateNonHermitianHamiltonian();
}


Mat2 RabiMCWF::GenerateHamiltonian(){
    return aDagMinusA * (this->eta * I);
}
Mat2 RabiMCWF::GenerateNonHermitianHamiltonian(){
    Mat2 jump_part(0,0,0,0);

    for(int i=0; i<2; i++){
        Mat2 jump_H = (this->jumpOperators[i]).Hermitian();
        jump_part += DOT(jump_H, this->jumpOperators[i]);
    }
    return this->ham - jump_part * ( 0.5 * this->hbar * I);
}
std::vector<Mat2> RabiMCWF::GenerateJumpOperators(){

    std::vector<Mat2> jumps(2);

    Mat2 j_up = a * sqrt(2 * this->kappa * (this->n_th));
    Mat2 j_down = aDag * sqrt(2 * this->kappa * (this->n_th + 1));

    jumps[0] = j_up;
    jumps[1] = j_down;
    return jumps;
}

void RabiMCWF::NoJump(double dt, double dp, int it, std::vector<Vec2>& phi){
    Mat2 lhs = Id - (this->hamNH * (I * dt/this->hbar));
    Vec2 res = lhs.dot(phi[it-1]);
    res /= sqrt(1-dp);
    phi[it] = res;
}
void RabiMCWF::Jump(double dt, double p_k, Mat2& jumpOp, int it, std::vector<Vec2>& phi){
    double mult = sqrt(dt/p_k);
    Vec2 rhs = (jumpOp.dot(phi[it-1]));
    phi[it] = rhs * mult;
}

Vec2R RabiMCWF::CalcPk(Vec2& phi){
    

    double y  = real(phi.Conj().dot((this->jumpOperators[0]).Hermitian().dot(this->jumpOperators[0].dot(phi))));
    double x  = real(phi.Conj().dot((this->jumpOperators[1]).Hermitian().dot(this->jumpOperators[1].dot(phi))));

    return Vec2R(x,y);
}



void RabiMCWF::Propagate(double dt, int nSteps, std::vector<Vec2>& phi, int measureDelta, int measureLength){
    phi[0].x = 1;
    phi[0].y = 0;
    int last_delta = 0;

    for(int it=1; it<nSteps; it++){

        
        Vec2R pK_all = this->CalcPk(phi[it-1]) * dt;
        double dp = pK_all.x + pK_all.y;
        double epsilon = MyRand();

        bool measure=false;

        if(it - last_delta > measureDelta){
            measure = true;
            if(it - measureDelta > measureDelta + measureLength){
                last_delta = it;
                measure = false;
            }
        }

        if(measure){
            double p_current = 0;
            double r = MyRand() * dp;

            bool fin= false;
            for(int i=0; i<2; i++){
                if(fin){ continue; }

                double pk = pK_all[i];
                p_current += pk;
                if(p_current > r){
                    this->Jump(dt, pk, this->jumpOperators[i], it, phi);
                    phi[it].normalize();
                    fin=true;
                }
            }

        }else{
            this->NoJump(dt, dp, it, phi);
            phi[it].normalize();
        }

    }

}


std::vector<Vec2R> RabiMCWF::RunExperiment(int nSteps, int measureDelta, int measureLength, double dt, int iterations){
    //Create vector to store our results
    std::vector<Vec2R> results(nSteps);

    std::vector<Vec2> phi(nSteps);

    for(int is=0; is<iterations; is++){
        //calculate all our phi
        this->Propagate(dt, nSteps, phi, measureDelta, measureLength);

        for(int it=0; it<nSteps; it++){
            results[it].x += real(conj(phi[it].x) * phi[it].x) / (1.0*iterations);
            results[it].y += real(conj(phi[it].y) * phi[it].y) / (1.0*iterations);
        }

    }

    return results;
}