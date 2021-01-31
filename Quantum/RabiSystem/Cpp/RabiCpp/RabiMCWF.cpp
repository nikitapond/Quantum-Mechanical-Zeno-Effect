#pragma once

#include "RabiMCWF.h"
#include "Vec2.h"
#include "Mat2.h"
#include <random>
double MyRand(){
    return ((double)rand())/(RAND_MAX+1);
}




Mat2 Sp(0, 0, 1, 0); //Sigma+
Mat2 Sm(0, 1, 0, 0);
Mat2 Ee(0, 0, 0, 1);
Mat2 Ve(0, 1, 1, 0);
Mat2 Id(1, 0, 0, 1);

Mat2 aDag(0, 0, 1, 0);
Mat2 a(0, 1, 0, 0);

Mat2 aDagMinusA(0, -1, 1, 0);
Mat2 aDagMultA(0, 0, 0, 1);

RabiMCWF::RabiMCWF() {
    hbar = 1.;
    n_th = 1.;
    kappa = 1.;
    eta = 1.;
    Init();
}
RabiMCWF::RabiMCWF(double hbar, double n_th, double kappa, double eta){
    this->hbar = hbar;
    this->n_th = n_th;
    this->kappa = kappa;
    this->eta = eta;
    Init();
}

void RabiMCWF::Init() {
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

        //std::cout << "herm: \n" << jump_H.ToString() << std::endl;
        //std::cout << "norm \n" << jumpOperators[i].ToString() << std::endl;

       // Mat2 temp = DOT(jump_H, this->jumpOperators[i]);
        jump_part += DOT(jump_H, this->jumpOperators[i]);

    //    std::cout << "jpart" << i << "\n" << temp.ToString() << std::endl;
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
   // std::cout << "lhs: \n" << lhs.ToString() << std::endl;
    //std::cout << "Matrix: \n" << lhs.ToString() << std::endl;
    //std::cout << "dp: " << dp << std::endl;
    Vec2 res = lhs.dot(phi[it-1]);
    res /= sqrt(1-dp);
    //std::cout << "Phi no jump: " << res.ToString() << std::endl;
    phi[it] = res;
}
void RabiMCWF::Jump(double dt, double p_k, Mat2& jumpOp, int it, std::vector<Vec2>& phi){
    double mult = sqrt(dt/p_k);
    //std::cout << "mult: " << mult << std::endl;

   // std::cout << "jump: \n" << jumpOp.ToString() << std::endl;

    Vec2 rhs = (jumpOp.dot(phi[it-1]));
    //std::cout << "Phi:\n" << phi[it - 1].ToString() << std::endl;
    //std::cout << "rhs\n" << rhs.ToString() << std::endl;
    phi[it] = rhs * mult;
   // std::cout << "fin\n" << phi[it].ToString() << std::endl;
}

Vec2R RabiMCWF::CalcPk(Vec2& phi){
    


   // Vec2 rhs_1 = phi.Conj().dot(this->jumpOperators[0].Hermitian());
   // Vec2 rhs_2 = rhs_1.dot(this->jumpOperators[0]);
    //double fin_y = real(rhs_2.dot(phi));

    double y = real(phi.Conj().dot(this->jumpOperators[0].Hermitian()).dot(this->jumpOperators[0]).dot(phi));
    double x = real(phi.Conj().dot(this->jumpOperators[1].Hermitian()).dot(this->jumpOperators[1]).dot(phi));



   // Vec2 rhs_1_x = phi.Conj().dot(this->jumpOperators[1].Hermitian());
    //Vec2 rhs_2_x = rhs_1_x.dot(this->jumpOperators[1]);
    //double fin_x = real(rhs_2_x.dot(phi));
    //std::cout << "fin x: " << fin_x << "\n" << std::endl;

    //double x  = real(phi.Conj().dot((this->jumpOperators[1]).Hermitian().dot(this->jumpOperators[1].dot(phi))));
   // double x = fin_x;
   // if (isnan(x) || isnan(y)) {
   //     std::cout << "yeah this is the issue\n x: " << x << "\t, y: " << y << std::endl;
    //}

    //std::cout << "PK shit: " << std::endl;
    //std::cout << "PK : " << x << "," << y << std::endl;
    //std::cout << "modified y " << fin_y << std::endl;
    //std::cout << phi.ToString() << std::endl;

    //std::cout << "rhs 1:\t" << rhs_y_1.ToString() << std::endl;
    //std::cout << "rhs 2:\t" << rhs_y_2.ToString() << std::endl;
    //std::cout << "y fin:\t" << y_fin << std::endl;
    return Vec2R(x,y);
}



bool RabiMCWF::Propagate(double dt, int nSteps, std::vector<Vec2>& phi, int measureDelta, int measureLength){
    phi[0].x = 1;
    phi[0].y = 0;
    int last_delta = 0;

    for(int it=1; it<nSteps; it++){
   
        Vec2R pK_all = this->CalcPk(phi[it-1]) * dt;


        double dp = pK_all.x + pK_all.y;


            

        double epsilon = MyRand();
        
        bool measure=false;

        if(it - last_delta >= measureDelta){
            measure = true;
            /*if(it - measureDelta > measureDelta + measureLength){
                last_delta = it;
                measure = false;
            }*/
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

            if (!fin) {

                std::cout << "Jump has failed!\n  r: " << r << "\nPk_all: " << pK_all.ToString()  << "\n r: " << r << std::endl;
                this->NoJump(dt, dp, it, phi);
                //  std::cout << "Phi: " << phi[it].x << "," << phi[it].y << std::endl;

                phi[it].normalize();
                return false;
                //std::cout << "Jump failed!, our result was " << phi[it].ToString() << std::endl;
            }
            measure = false;

        }else{
           // std::cout << "no measure" << std::endl;
            this->NoJump(dt, dp, it, phi);
          //  std::cout << "Phi: " << phi[it].x << "," << phi[it].y << std::endl;

            phi[it].normalize();
            //std::cout << "After no jump, our result was " << phi[it].ToString() << std::endl;

          //  std::cout << "Phi: " << phi[it].x << "," << phi[it].y << std::endl;
        }

    }

}


std::vector<Vec2R> RabiMCWF::RunExperiment(int nSteps, int measureDelta, int measureLength, double dt, int iterations){

    /*std::cout <<" J0:\n" << this->jumpOperators[0].ToString() << std::endl;
    std::cout <<" J1:\n" <<  this->jumpOperators[1].ToString() << std::endl;

    std::cout << "ham\n"<< this->ham.ToString() << std::endl;
    std::cout << "hamNH\n" << this->hamNH.ToString() << std::endl;
    *///Create vector to store our results
    std::vector<Vec2R> results(nSteps);

    std::vector<Vec2> phi(nSteps);
    

    srand(time(NULL));

    float progress = 0;
    for(int is=0; is<iterations; is++){
        //srand(time(NULL));

        //calculate all our phi
        if (!this->Propagate(dt, nSteps, phi, measureDelta, measureLength))
            return results;
        
        for(int it=0; it<nSteps; it++){
            results[it].x += (real(conj(phi[it].x) * phi[it].x));
            results[it].y += (real(conj(phi[it].y) * phi[it].y));
        }
        progress = (is * 1.0) / iterations;
        std::cout << "Progress: " << int(progress * 100) << "%\r";
        std::cout.flush();
    }
    for (int it = 0; it < nSteps; it++) {
        results[it].x /= (1.0 * iterations);
        results[it].y /= (1.0 * iterations);
    }
    
    return results;
}