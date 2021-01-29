

#include "RabiMCWF.h"

// #include <iostream>
// #include <complex>
// #include <cstdlib>
// #include <cstdio>
// #include <vector>
// const int N = 2;

// const double Sp[N][N] = { {0., 0.},
//                           {1., 0.} }; // Sigma+
// const double Sm[N][N] = { {0., 1.},
//                           {0., 0.} }; // Sigma-
// const double Ee[N][N] = { {0., 0.},
//                           {0., 1.} }; // Sigma+ * Sigma-
// const double Ve[N][N] = { {0., 1.},
//                           {1., 0.} }; // Sigma+ + Sigma-
// const double Id[N][N] = { {1., 0.},
//                           {0., 1.} }; // Identity

// #define complex std::complex<double>
// const complex I = complex(0., 1.);

// class Mat2{

// private:
//     complex m11, m12, m21, m22;
// public:
//     Mat2(){
//         m11=m12=m21=m22= I * 0.0;
//     }
//     Mat2(complex m11, complex m12, complex m21, complex m22){
//         this->m11 = m11;
//         this->m12 = m12;
//         this->m21 = m21;
//         this->m22 = m22;
//     }

//     Mat2 dot(Mat2 &b){
//         complex m11 = m11 * b.m11 + m12 * b.m21;
//         complex m12 = m11 * b.m12 + m12 * b.m22;

//         complex m21 = m21 * b.m11 + m22 * b.m21;
//         complex m22 = m21 * b.m12 + m22 * b.m22;

//         return Mat2(m11, m12, m21, m22);
//     };
// };


// inline Mat2 DOT(Mat2 &a, Mat2 &b){
//     return a.dot(b);
// };

// class RabiMCWF{

// private:
//     float hbar, n_th, kappa, eta;

//     Mat2 ham; //Hamiltonian for this system
//     Mat2 hamNH; //Non hermitian hamiltonian
//     std::vector<Mat2> jumpOperators; //Jump operators

//     Mat2 GenerateHamiltonian();
//     Mat2 GenerateNonHermitianHamiltonian();
//     std::vector<Mat2> GenerateJumpOperators();
// public:
//     RabiMCWF(){
//         hbar = 1.;
//         n_th = 1.;
//         kappa = 1.;
//         eta = 1.;
//     }
//     RabiMCWF(float hbar, float n_th, float kappa, float eta);

//     std::vector<float> RunExperiment(int nSteps, float measureDelta, float measureLength, float dt, int iterations);

// };

RabiMCWF::RabiMCWF(float hbar, float n_th, float kappa, float eta){
    this->hbar = hbar;
    this->n_th = n_th;
    this->kappa = kappa;
    this->eta = eta;
}

std::vector<float> RabiMCWF::RunExperiment(int nSteps, float measureDelta, float measureLength, float dt, int iterations){
    std::vector<float> results(nSteps);


    return results;
}