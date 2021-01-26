const double Omega = 0.06;
const double delta = 0.01;
const double Gamma = 0.02;
const long nt = 4000;
const long ns = 100;
const double dt = 0.1;

const int N = 2;

const double Sp[N][N] = { {0., 0.},
                          {1., 0.} }; // Sigma+
const double Sm[N][N] = { {0., 1.},
                          {0., 0.} }; // Sigma-
const double Ee[N][N] = { {0., 0.},
                          {0., 1.} }; // Sigma+ * Sigma-
const double Ve[N][N] = { {0., 1.},
                          {1., 0.} }; // Sigma+ + Sigma-
const double Id[N][N] = { {1., 0.},
                          {0., 1.} }; // Identity

#include <complex>
#define complex std::complex<double>
const complex I = complex(0., 1.);

#include <cstdlib>
#include <cstdio>

void generateHamiltonian(complex hamil[N][N]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            hamil[i][j] = -delta * Ee[i][j] + 0.5 * Omega * Ve[i][j] - 0.5 * I * Gamma * Ee[i][j];
}

void generateEvolutionU(complex hamil[N][N], double dp, complex U[N][N]) {
    double mu = 1. / sqrt(1. - dp);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            U[i][j] = mu * (Id[i][j] - I * dt * hamil[i][j]);
}

void normalize(complex phi[N]) {
    double norm = 0.;
    for (int i = 0; i < N; i++)
        norm += real(conj(phi[i]) * phi[i]);
    norm = sqrt(norm);
    for (int i = 0; i < N; i++)
        phi[i] /= norm;
}

void propagate(complex U[N][N], complex* phi_old, complex* phi_new) {
    for (int i = 0; i < N; i++) {
        phi_new[i] = 0.;
        for (int k = 0; k < N; k++)
            phi_new[i] += U[i][k] * phi_old[k];
    }
}

void propagate_single_traj(complex hamil[N][N], complex phi[nt][N]) {
    phi[0][0] = 1.; phi[0][1] = 0.;
    for (long it = 1; it < nt; it++) {
        double dp = real(Gamma * dt * conj(phi[it - 1][1]) * phi[it - 1][1]); // Gamma*dt*<phi(t)|S^+ S^-|phi(t)>
        double epsilon = (double)rand() / RAND_MAX;
        if (epsilon >= dp) {
            complex U[N][N];
            generateEvolutionU(hamil, dp, U);
            propagate(U, phi[it - 1], phi[it]);
            normalize(phi[it]);
        }
        else {
            phi[it][0] = 1.;
            phi[it][1] = 0.;
        }
    }
}

void dump_population_to_file(double population[nt][N]) {
    FILE* file = fopen("C:/Users/nikit/github/ZenoQM/Zeno_QM/Python/population.dat", "w");
    for (long it = 0; it < nt; it++)
        fprintf(file, "%le %le %le\n", (0. + it * dt) /** Gamma*/, population[it][0], population[it][1]);
    fclose(file);
}

void dump_phi_to_file(complex phi[nt][N]) {
    FILE* file = fopen("wf.dat", "w");
    for (long it = 0; it < nt; it++)
        fprintf(file, "%le %le %le %le %le\n", (0. + it * dt) /** Gamma*/,
            real(phi[it][0]), real(phi[it][0]), real(phi[it][1]), real(phi[it][1]));
    fclose(file);
}

int main3(int argc, char* argv[]) {
    complex hamil[N][N];
    generateHamiltonian(hamil);

    complex phi[nt][N];
    double population[nt][N];
    for (long it = 0; it < nt; it++)
        for (int i = 0; i < N; i++)
            population[it][i] = 0.;

    for (int is = 0; is < ns; is++) {
        srand(is);
        propagate_single_traj(hamil, phi);
        for (long it = 0; it < nt; it++)
            for (int i = 0; i < N; i++)
                population[it][i] += real(conj(phi[it][i]) * phi[it][i]) / (1. * ns);
    }

    dump_population_to_file(population);
    return 0;
}