#include "HfSolver.hpp"
#include "Molecule.hpp"
#include "Helper.hpp"
#include <fstream>
#include <algorithm>

#define MAXITER 100
#define DELTA_1 1e-12
#define DELTA_2 1e-11

#define INDEX(i, j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

HfSolver::HfSolver(Molecule &m, bool toprint):
    EnergySolver(m){
    this -> toprint = toprint;
    computed = false;
    initialize();
}

void HfSolver::read_dipole(double ****mux, double ****muy, double ****muz){
    std::string path = dir + "/mux.dat";
    mux = readMatrix(path.c_str());

    path = dir + "/muy.dat";
    muy = readMatrix(path.c_str());

    path = dir + "/muz.dat";
    muz = readMatrix(path.c_str());
}

void HfSolver::compute_dipole(){
    double ****mux, ****muy, ****muz;
    read_dipole(mux, muy, muz);

    double dx = 0.0, dy = 0.0, dz = 0.0;
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            dx += 2 * D(i, j) * mux[i][j];
            dy += 2 * D(i, j) * muy[i][j];
            dz += 2 * D(i, j) * muz[i][j];
        }
    }

    for (int i = 0; i < natom; i++){
        dx += atoms[i].zval * atoms[i].x;
        dy += atoms[i].zval * atoms[i].y;
        dz += atoms[i].zval * atoms[i].z;
    }
    printf("Mux = %20.12f\nMuy = %20.12f\nMuz = %20.12f\n", dx, dy, dz);
    printf("Total dipole moment (au) = %20.12f\n", dx + dy + dz);
}

double HfSolver::calc_hf_energy(){
    double E = 0.0;
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            E += D(i, j) * (ham[i][j] + F(i, j));
        }
    }
    return E;
}

double HfSolver::compute(){

    if (computed && toprint) initialize(); 

    double E_prev = 0.0;

    double E_curr = calc_hf_energy();

    double delta_E = E_curr - E_prev;
    double rms = 1.0;

    int count = 0;

    if (toprint){ 
        printf("Iter\t\tE(elec)  \t\tE(tot)  \t\tDelta(E)  \t\tRMS(D)\n");
        printf("%02d%21.12f%21.12f\n", count, E_curr, E_curr + enuc);
    }

    while (count < MAXITER && (abs(delta_E) >= DELTA_1 || rms >= DELTA_2)){
        E_prev = E_curr;

        updateFock();

        Matrix new_D = Matrix::Zero(norb, norb);
        updateDensity(new_D);

        if(toprint && count == 0){
            printf("\tFock Matrix:\n\n");
            Helper::print_matrix(F);
        }

        rms = Helper::calc_rms(D, new_D);
        
        D = new_D;
        E_curr = calc_hf_energy();
        delta_E = E_curr - E_prev;

        count++;
        if(toprint){
            printf("%02d%21.12f%21.12f%21.12f%21.12f\n", count, E_curr, 
                                                     E_curr + enuc, delta_E, rms); 
        }
    }

    if (toprint){
        Helper::print_matrix(D);
        compute_dipole();
    }

    computed = true;

    return E_curr;
}

void HfSolver::initialize(){

    // initialize S^(-1/2) matrix
    Matrix S = Matrix::Zero(norb, norb);

    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            S(i, j) = s[i][j];
        }
    }

    Eigen::SelfAdjointEigenSolver<Matrix> solver(norb);

    solver.compute(S);
    Matrix evecs_S = solver.eigenvectors();
    Matrix evals_S = solver.eigenvalues();

    Matrix A = evals_S.asDiagonal();

    Matrix isqrt_A = A.inverse().cwiseSqrt();

    isqrt_S = evecs_S * isqrt_A * evecs_S.transpose();

    F = Matrix::Zero(norb, norb);
    // initialize F matrix
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            F(i, j) = ham[i][j];
        }
    }

    // initialize C matrix
    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;
    
    solver.compute(Fp);
    Matrix Cp = solver.eigenvectors();

    C = isqrt_S * Cp;

    D = Matrix::Zero(norb, norb);
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            double sum = 0.0;
            for (int k = 0; k < nomo; k++){
                sum += C(i, k) * C(j, k);
            }
            D(i, j) = sum;
        }
    }

    if (toprint){
        printf("\tS^-1/2 Matrix:\n\n");
        Helper::print_matrix(isqrt_S);

        printf("\tInitial F' Matrix:\n\n");
        Helper::print_matrix(Fp);

        printf("\tInitial C Matrix:\n\n");
        Helper::print_matrix(C);
        
        printf("\tInitial Density Matrix:\n\n");
        Helper::print_matrix(D);
    }
}

void HfSolver::updateFock(){
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            
            F(i, j) = ham[i][j];
            
            for (int k = 0; k < norb; k++){
                for (int l = 0; l < norb; l++){
                    int ij = INDEX(i, j);
                    int kl = INDEX(k, l);
                    int ik = INDEX(i, k);
                    int jl = INDEX(j, l);

                    F(i, j) += D(k, l) * (2 * eri[INDEX(ij, kl)] - eri[INDEX(ik, jl)]);
                }
            }
        }
    }
}

void HfSolver::updateDensity(Matrix &new_D){
    
    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(Fp);
    Matrix evecs_Fp = solver.eigenvectors();
    Matrix evals_Fp = solver.eigenvalues();

    C = isqrt_S * evecs_Fp;

    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            double sum = 0.0;
            for (int k = 0; k < nomo; k++){
                sum += C(i, k) * C(j, k);
            }
            new_D(i, j) = sum;
        }
    }
}