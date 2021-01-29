#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>
const int N = 2;

const Mat2 Sp(0,0,1,0); //Sigma+
const Mat2 Sm(0,1,0,0);
const Mat2 Ee(0,0,1,0);
const Mat2 Ve(0,1,1,0);
const Mat2 Id(1,0,1,0);

const Mat2 aDag(0,0,1,0);
const Mat2 a(0,1,0,0);

const Mat2 aDagMinusA(0,-1,1,0);
const Mat2 aDagMultA(0,0,0,1);

#define complex std::complex<double>
const complex I = complex(0., 1.);

class Mat2{

private:
    complex m11, m12, m21, m22;
public:
    Mat2(){
        m11=m12=m21=m22= I * 0.0;
    }
    Mat2(complex m11, complex m12, complex m21, complex m22){
        this->m11 = m11;
        this->m12 = m12;
        this->m21 = m21;
        this->m22 = m22;
    }


    Mat2 Hermitian(){
        complex m11 = conj(this->m11);
        complex m12 = conj(this->m21);
        complex m21 = conj(this->m12);
        complex m22 = conj(this->m22);

        return Mat2(m11, m12, m21, m22);
    }

    Mat2 dot(Mat2 &b){
        complex m11 = m11 * b.m11 + m12 * b.m21;
        complex m12 = m11 * b.m12 + m12 * b.m22;

        complex m21 = m21 * b.m11 + m22 * b.m21;
        complex m22 = m21 * b.m12 + m22 * b.m22;

        return Mat2(m11, m12, m21, m22);
    };
    Mat2 operator*(complex c){
        return Mat2(this->m11 * c, this->m12 * c, this->m21 * c, this->m22 * c);
    }

    Mat2& operator*=(complex c){
        this->m11 = this->m11 * c;
        this->m12 = this->m12 * c;
        this->m21 = this->m21 * c;
        this->m22 = this->m22 * c;
        return *this;
    }

    Mat2 operator*(double c){
        return Mat2(this->m11 * c, this->m12 * c, this->m21 * c, this->m22 * c);
    }
    Mat2& operator*=(double c){
        this->m11 = this->m11 * c;
        this->m12 = this->m12 * c;
        this->m21 = this->m21 * c;
        this->m22 = this->m22 * c;
        return *this;
        
    }

    Mat2 operator+(Mat2 c){
        complex m11 = this->m11 + c.m11;
        complex m12 = this->m12 + c.m12;
        complex m21 = this->m21 + c.m21;
        complex m22 = this->m22 + c.m22;
        return Mat2(m11,m12,m21,m22);
    }

    Mat2& operator+=(Mat2 c){
        this->m11 = this->m11 + c.m11;
        this->m12 = this->m12 + c.m12;
        this->m21 = this->m21 + c.m21;
        this->m22 = this->m22 + c.m22;
        return *this;
    }
    Mat2 operator-(Mat2 c){
        complex m11 = this->m11 - c.m11;
        complex m12 = this->m12 - c.m12;
        complex m21 = this->m21 - c.m21;
        complex m22 = this->m22 - c.m22;
        return Mat2(m11,m12,m21,m22);
    }
};

inline Mat2 DOT(Mat2 &a, Mat2 &b){
    return a.dot(b);
};


class RabiMCWF{

private:
    float hbar, n_th, kappa, eta;

    Mat2 ham; //Hamiltonian for this system
    Mat2 hamNH; //Non hermitian hamiltonian
    std::vector<Mat2> jumpOperators; //Jump operators

    Mat2 GenerateHamiltonian();
    Mat2 GenerateNonHermitianHamiltonian();
    std::vector<Mat2> GenerateJumpOperators();

    void Normalize(complex phi[2]);
    void NoJump(float dt, float dp, int it, complex phi[2]);
public:
    RabiMCWF(){
        hbar = 1.;
        n_th = 1.;
        kappa = 1.;
        eta = 1.;
    }
    RabiMCWF(float hbar, float n_th, float kappa, float eta);

    std::vector<float> RunExperiment(int nSteps, float measureDelta, float measureLength, float dt, int iterations);

};



RabiMCWF::RabiMCWF(float hbar, float n_th, float kappa, float eta){
    this->hbar = hbar;
    this->n_th = n_th;
    this->kappa = kappa;
    this->eta = eta;
    this->ham = GenerateHamiltonian();
    this->jumpOperators = GenerateJumpOperators();
    this->hamNH = GenerateNonHermitianHamiltonian();
}


Mat2 RabiMCWF::GenerateHamiltonian(){
    return aDagMinusA * this->eta * I;
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

void RabiMCWF::Normalize(complex phi[2]){
    double norm = 0;
    for(int i=0; i<2; i++){
        norm += real(conj(phi[i]) * phi[i]);
    }
    for(int i=0; i<2; i++){
        phi[i] /= norm;
    }
}
void RabiMCWF::NoJump(float dt, float dp, int it, )



std::vector<float> RabiMCWF::RunExperiment(int nSteps, float measureDelta, float measureLength, float dt, int iterations){
    std::vector<float> results(nSteps);


    return results;
}