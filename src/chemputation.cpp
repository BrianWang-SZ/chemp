#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include "molecule.hpp"
#include "masses.h"
#include "assert.h"

#define MAXORB 100
#define MAXITER 100
#define MAXERR 8
#define DELTA_1 1e-12
#define DELTA_2 1e-11

#define INDEX(i, j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))
#define here printf("here");

Molecule::~Molecule(){
    delete[] atoms;
    delete[] eri;
    delete[] ioff;

    free2dvec(s, norb);
    free2dvec(v, norb);
    free2dvec(t, norb);

    free2dvec(ham, norb);
    free2dvec(mux, norb);
    free2dvec(muy, norb);
    free2dvec(muz, norb);
}

// lookup array

void Molecule::lookupTable(){
    int* ioff = new int[MAXORB];

    ioff[0] = 0;

    for (int i = 1; i < MAXORB; i++){
        ioff[i] = ioff[i - 1] + i;
    }

    this -> ioff = ioff;
}


void Molecule::hamiltonian(){
    double **ham = create2dvec(norb);
    
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            ham[i][j] = t[i][j] + v[i][j];
        }
    }
    this -> ham = ham;

    printf("\tCore Hamiltonian:\n\n");
    print_matrix(ham, norb);

}

void Molecule::_one_electron(const char *dir){

    FILE *in;
    std::string direc = dir;
    std::string path = direc + "/enuc.dat";
    
    if ((in = fopen(path.c_str(), "r")) == NULL) {
        perror("Error opening enuc.dat\n");
        exit(-1);
    }

    double enuc;
    if (fscanf(in, " %lf", &enuc) != 1) {
        perror("Error reading nuclear energy\n");
        exit(-1);
    }
    this -> enuc = enuc;

    if (fclose(in)){
        perror("Error closing enuc.dat\n");
        exit(-1);
    }

    printf("Nuclear repulsion energy = %20.10f\n", enuc);

    path = direc + "/s.dat";

    this -> s = readMatrix(path.c_str());
    printf("\n\tOverlap Integrals:\n\n");
    print_matrix(s, norb);
    
    path = direc + "/t.dat";
    this -> t = readMatrix(path.c_str());
    printf("\tKinetic-Energy Integrals:\n\n");
    print_matrix(t, norb);

    path = direc + "/v.dat";
    this -> v = readMatrix(path.c_str());
    printf("\tNuclear Attraction Integrals:\n\n");
    print_matrix(v, norb);
}

void Molecule::read_two_electron(const char *dir){
    FILE *in;
    std::string direc = dir;
    std::string path = direc + "/eri.dat";
    
    if ((in = fopen(path.c_str(), "r")) == NULL) {
        perror("Error opening eri.dat");
        exit(-1);
    }

    int i = norb - 1, j = norb - 1, 
        k = norb - 1, l = norb - 1;

    double val = 0.0;

    int ij = INDEX(i, j);
    int kl = INDEX(k, l);
    int ijkl = INDEX(ij, kl);

    double *eri = new double[ijkl];

    for(int i = 0; i < ijkl; i++){
        eri[i] = 0.0;
    }

    while (fscanf(in, " %d %d %d %d %lf", &i, &j, &k, &l, &val) == 5) {
        i--; j--; k--; l--;
        ij = INDEX(i, j);
        kl = INDEX(k, l);
        ijkl = INDEX(ij, kl);
        eri[ijkl] = val;
    }

    if (fclose(in)){
        perror("Error closing eri.dat\n");
        exit(-1);
    }
    this -> eri = eri;
}

double** Molecule::readMatrix(const char* path){

    if (this -> norb == 0) {
        int norb = calc_norb(path);
        if (norb == -1){
            perror("Exceeds maximum orbitals!");
            exit(-1);
        }
        this -> norb = norb;
    }
    
    FILE *in;
    double **mat = create2dvec(norb);

    if ((in = fopen(path, "r")) == NULL) {
        fprintf(stderr, "Error opening %s\n", path);
        perror("");
        exit(-1);
    }
    
    int a, b;
    double c;
    int count = 0;

    while (fscanf(in, "%d %d %lf", &a, &b, &c) == 3) {
        a--; b--;
        mat[a][b] = c;
        if (a != b) mat[b][a] = mat[a][b];
        count++;
    }
    
    if (calc_norb(path) != norb) {
        fprintf(stderr, "Inconsistent input file %s\n", path);
        perror("");
        exit(-1);
    }

    if (fclose(in)){
        fprintf(stderr, "Error closing %s\n", path);
        perror("");
        exit(-1);
    }
    return mat;
}

Molecule::Molecule(const char *dir){
    lookupTable();
    read_coord(dir);
    read_one_electron(dir);
    hamiltonian();
    read_two_electron(dir);
    read_dipole(dir);
    compute();
}

void Molecule::print_matrix(double **mat, int size){

    for (int i = 0; i < size; i++){
        printf("%12d", i + 1);
    }

    printf("\n\n");

    for (int i = 0; i < size; i++){
        printf("%5d", i + 1);
        for (int j = 0; j < size; j++){
            printf("%12.7f", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void Molecule::print_matrix(Matrix mat){
    for (int i = 0; i < mat.cols(); i++){
        printf("%12d", i + 1);
    }

    printf("\n\n");

    for (int i = 0; i < mat.rows(); i++){
        printf("%5d", i + 1);
        for (int j = 0; j < mat.cols(); j++){
            printf("%12.7f", mat(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

int Molecule::calc_norb(const char* path){
    // count number of lines in the file
    std::ifstream inFile(path); 
    int nline = std::count(std::istreambuf_iterator<char>(inFile), 
                            std::istreambuf_iterator<char>(), '\n');
    
    // calculate orbital number
    for (int i = 0; i <= MAXORB; i++){
        if (i + pow(i, 2) == nline * 2) {
            return i;
        }
    }
    return -1;
}

void Molecule::compute(){
    Matrix isqrt_S = Matrix::Zero(norb, norb);
    Matrix C = Matrix::Zero(norb, norb);
    Matrix D = Matrix::Zero(norb, norb);
    Matrix F = Matrix::Zero(norb, norb);

    initialize(0, isqrt_S, C, D, F);
    double E_scf = compute_hf(isqrt_S, C, D, F);
    double E_mp2 = compute_mp2(C, F, isqrt_S);
    printf("Escf = %20.12f\nEmp2 = %20.12f\nEtot = %20.12f\n", E_scf + enuc, 
                                                             E_mp2, E_scf + enuc + E_mp2);
    double E_cc = 0.0;
    
    if (( E_cc = compute_ccsd(C, F, isqrt_S)) != 0.0){
        printf("Etot = %21.12f\n", E_scf + E_cc + enuc);
    }
}

void Molecule::initialize(int toprint, Matrix &isqrt_S, Matrix &C, Matrix &D, Matrix &F){

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

    if (toprint){
        printf("\tS^-1/2 Matrix:\n\n");
        print_matrix(isqrt_S);
    }

    // initialize F matrix
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            F(i, j) = ham[i][j];
        }
    }

    // initialize C matrix
    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;
    
    if (toprint){
        printf("\tInitial F' Matrix:\n\n");
        print_matrix(Fp);
    }

    solver.compute(Fp);
    Matrix Cp = solver.eigenvectors();

    C = isqrt_S * Cp;

    if (toprint){
        printf("\tInitial C Matrix:\n\n");
        print_matrix(C);
    }
    
    // calculate number of occupied orbitals
    int nomo = calc_nomo();

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
        printf("\tInitial Density Matrix:\n\n");
        print_matrix(D);
    }
}

double Molecule::compute_hf(Matrix isqrt_S, Matrix &C, Matrix &D, Matrix &F){

    double E_prev = 0.0;

    double E_curr = calc_hf_energy(D, F);

    double delta_E = E_curr - E_prev;
    double rms = 1.0;

    int count = 0;

    printf("Iter\t\tE(elec)  \t\tE(tot)  \t\tDelta(E)  \t\tRMS(D)\n");
    printf("%02d%21.12f%21.12f\n", count, E_curr, E_curr + enuc);

    while (count < MAXITER && (abs(delta_E) > DELTA_1 || rms > DELTA_2)){
        E_prev = E_curr;

        updateFock(F, D);

        Matrix new_D = Matrix::Zero(norb, norb);
        updateDensity(new_D, isqrt_S, C, F);

        if(count == 0){
            printf("\tFock Matrix:\n\n");
            print_matrix(F);
        }

        rms = calc_rms(D, new_D);
        
        D = new_D;
        E_curr = calc_hf_energy(D, F);
        delta_E = E_curr - E_prev;

        count++;
        printf("%02d%21.12f%21.12f%21.12f%21.12f\n", count, E_curr, 
                                                     E_curr + enuc, delta_E, rms); 
    }

    print_matrix(D);
    return E_curr;
    //compute_dipole(D);
}

double compute_diis(Matrix D, Matrix F, Matrix isqrt_S){
    //Matrix *err = new Matrix[MAXERR];
    //Matrix *f = new Matrix[MAXERR];

    //Matrix e = isqrt_S.transpose() * (F * D * S - S * D * F) * isqrt_S;
        
        // if (count < MAXERR){
        //     err[count] = e;
        //     f[count] = F;
        // } else {
        //     for (int i = 1; i < MAXERR; i++){
        //         err[i - 1] = err[i];
        //         f[i - 1] = f[i];
        //     }

        //     err[MAXERR - 1] = e;
        //     f[MAXERR - 1] = F;
        // }
        
        //compute new density matrix
        //int max = (count > MAXERR) ? MAXERR : count;
        
        // if (count > 2){

            // Eigen::MatrixXd B(max + 1, max + 1);

            // for (int i = 0; i < max; i++){
            //     for (int j = 0; j < max; j++){
            //         B(i, j) = (err[i] * err[j].transpose()).trace();
            //     }
            // }

            // for (int i = 0; i < max; i++){
            //     B(max, i) = -1;
            //     B(i, max) = -1;
            // }

            // B(max, max) = 0;

            // Eigen::VectorXd b(max + 1);
            // b(max) = -1;
            
            // Eigen::VectorXd c = B.householderQr().solve(b);

            // Matrix Fp(norb, norb);

            // for (int i = 0; i < max; i++){
            //     Fp += c(i) * f[i];
            // }
            //updateDensity(new_D, isqrt_S, Fp, C);
            
        //} else {

        // }
        return 0.0;
}

double* Molecule::spatial_atom(double *eri, Matrix C){
    int max = INDEX(norb, norb);
    double *moeri = new double[INDEX(max, max)];

    for (int i = 0; i < INDEX(max, max); i++){
        moeri[i] = 0.0;
    }

    double ****M = create4dmat(norb);
    double ****N = create4dmat(norb);
    double ****P = create4dmat(norb);
    double ****Q = create4dmat(norb);
    
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            for (int k = 0; k < norb; k++){
                for(int l = 0; l < norb; l++){
                    for (int m = 0; m < norb; m++){
                        int ij = INDEX(i, j);
                        int kl = INDEX(k, l);
                        M[m][j][k][l] += C(i, m) * eri[INDEX(ij, kl)];
                    }
                }
            }
        }
    }

    for (int m = 0; m < norb; m++){
        for (int j = 0; j < norb; j++){
            for (int k = 0; k < norb; k++){
                for(int l = 0; l < norb; l++){
                    for (int n = 0; n < norb; n++){
                        N[m][n][k][l] += C(j, n) * M[m][j][k][l];
                    }
                }
            }
        }
    }

    for (int m = 0; m < norb; m++){
        for (int n = 0; n < norb; n++){
            for (int k = 0; k < norb; k++){
                for(int l = 0; l < norb; l++){
                    for (int p = 0; p < norb; p++){
                        P[m][n][p][l] += C(k, p) * N[m][n][k][l];
                    }
                }
            }
        }
    }

    for (int m = 0; m < norb; m++){
        for (int n = 0; n < norb; n++){
            for (int p = 0; p < norb; p++){
                for(int l = 0; l < norb; l++){
                    for (int q = 0; q < norb; q++){
                        Q[m][n][p][q] += C(l, q) * P[m][n][p][l];
                    }
                }
            }
        }
    }

    for (int m = 0; m < norb; m++){
        for (int n = 0; n < norb; n++){
            for (int p = 0; p < norb; p++){
                for (int q = 0; q < norb; q++){
                    int mn = INDEX(m, n);
                    int pq = INDEX(p, q);
                    int mnpq = INDEX(mn, pq);
                    moeri[mnpq] = Q[m][n][p][q];
                }
            }
        }
    }

    free4dmat(M, norb);
    free4dmat(N, norb);
    free4dmat(P, norb);
    free4dmat(Q, norb);

    return moeri;
}

int Molecule::calc_nomo(){
    int nelec = 0;

    for (int i = 0; i < natom; i++){
        nelec += atoms[i].zval;
    }

    return nelec / 2;
}

double Molecule::compute_mp2(Matrix C, Matrix F, Matrix isqrt_S){

    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(Fp);
    Matrix evals = solver.eigenvalues();

    double *moeri = spatial_atom(eri, C);

    int nomo = calc_nomo();

    double E_mp2 = 0.0;

    for (int i = 0; i < nomo; i++){
        for (int a = nomo; a < norb; a++){
            for (int j = 0; j < nomo; j++){
                for (int b = nomo; b < norb; b++){
                    int ia = INDEX(i, a);
                    int jb = INDEX(j, b);
                    int ib = INDEX(i, b);
                    int ja = INDEX(j, a);

                    E_mp2 += moeri[INDEX(ia, jb)] * (2 * moeri[INDEX(ia, jb)] - 
                                                        moeri[INDEX(ib, ja)]) /
                            (evals(i) + evals(j) - evals(a) - evals(b));
                }
            }
        }
    }
    return E_mp2;
}

double**** Molecule::create4dmat(int size){
    double ****mat = new double***[size];
    for (int i = 0; i < size; i++){
        mat[i] = new double**[size];
        for (int j = 0; j < size; j++){
            mat[i][j] = new double*[size];
            for (int k = 0; k < size; k++){
                mat[i][j][k] = new double[size];
                for (int l = 0; l < size; l++){
                    mat[i][j][k][l] = 0.0;
                }
            }
        } 
    }
    return mat;
}

void Molecule::free4dmat(double ****mat, int size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            for (int k = 0; k < size; k++){
                delete[] mat[i][j][k];
            }
            delete[] mat[i][j];
        }
        delete[] mat[i];
    }
    delete[] mat;
}

double Molecule::calc_rms(Matrix mat, Matrix newmat){
    if (mat.rows() != newmat.rows() || mat.cols() != newmat.cols()) return 0.0;
    
    double sum = 0.0;
    for (int i = 0; i < mat.rows(); i++){
        for (int j = 0; j < mat.cols(); j++){
            sum += pow(mat(i, j) - newmat(i, j), 2);
        }
    }
    return sqrt(sum);
}

double Molecule::calc_rms(double **mat, double **newmat, int size){
    
    double sum = 0.0;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            sum += pow(mat[i][j] - newmat[i][j], 2);
        }
    }
    return sqrt(sum);
}

double Molecule::calc_hf_energy(Matrix D, Matrix F){
    double E = 0.0;
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            E += D(i, j) * (ham[i][j] + F(i, j));
        }
    }
    return E;
}

void Molecule::updateFock(Matrix &F, Matrix D){
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

void Molecule::mobasis(Matrix C, Matrix F){
    Matrix moF = Matrix::Zero(norb, norb);
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            double sum = 0.0;
            for (int k = 0; k < norb; k++){
                for (int l = 0; l < norb; l++){
                    sum += C(k, j) * C(l, i) * F(k, l);
                }
            }
            F(i, j) = sum;
        }
    }
    print_matrix(moF);
}

void Molecule::updateDensity(Matrix &new_D, Matrix isqrt_S, Matrix &C, Matrix F){
    
    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(Fp);
    Matrix evecs_Fp = solver.eigenvectors();
    Matrix evals_Fp = solver.eigenvalues();

    C = isqrt_S * evecs_Fp;

    int nomo = calc_nomo();

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

void Molecule::read_coord(const char *dir){

    std::string direc = dir;
    std::string path = direc + "/geom.dat";
    
    // open input file
    FILE *in;
    
    if ((in = fopen(path.c_str(), "r")) == NULL) {
        perror("Error opening geom.dat");
        exit(-1);
    }
    
    // read in the number of atoms
    int num;
    if (fscanf(in, "%d", &num) != 1) {       
        perror("Error reading number of atoms");
        exit(-1);
    }
    natom = num;

    // read each atom
    int count = 0;
    
    Atom *atoms = new Atom[natom];

    double zval = 0.0, x = 0.0,
              y = 0.0, z = 0.0;
    
    while (fscanf(in, "%lf %lf %lf %lf", &zval, &x, &y, &z) == 4) {
        
        Atom *tmp = new Atom;
        tmp -> zval = (int) zval;
        tmp -> x = x;
        tmp -> y = y;
        tmp -> z = z;
        
        atoms[count] = *tmp; 
        count++;
    }

    if (count < num) {
        perror("Atom number does not match data");
        exit(-1);
    }

    if (fclose(in)){
        perror("Error closing file stream!");
        exit(-1);
    }

    this -> atoms = atoms;
}

void Molecule::read_dipole(const char *dir){
    std::string direc = dir;
    std::string path = direc + "/mux.dat";
    this -> mux = readMatrix(path.c_str());

    path = direc + "/muy.dat";
    this -> muy = readMatrix(path.c_str());

    path = direc + "/muz.dat";
    this -> muz = readMatrix(path.c_str());
}

double**** Molecule::spatial_to_spin(double *moeri){
    int max = 2 * norb;

    double ****mospin = create4dmat(max);

    for (int p = 0; p < max; p++){
        for (int q = 0; q < max; q++){
            for (int r = 0; r < max; r++){
                for (int s = 0; s < max; s++){
                    int pr = INDEX(p / 2, r / 2);
                    int qs = INDEX(q / 2, s / 2);
                    int is_same_parity = (p % 2 == r % 2) * (q % 2 == s % 2);
                    double a = moeri[INDEX(pr, qs)] * is_same_parity;
                    
                    int ps = INDEX(p / 2, s / 2);
                    int qr = INDEX(q / 2, r / 2);
                    is_same_parity = (p % 2 == s % 2) * (q % 2 == r % 2);
                    double b = moeri[INDEX(ps, qr)] * is_same_parity;
                    mospin[p][q][r][s] = a - b;
                }
            }
        }
    }
    return mospin;
}

double** Molecule::create2dvec(int size){
    double **vec = new double*[size];
    
    for (int i = 0; i < size; i++){
        vec[i] = new double[size];
        for (int j = 0; j < size; j++){
            vec[i][j] = 0.0;
        }
    }
    return vec;
}

void Molecule::free2dvec(double **vec, int size){
    for (int i = 0; i < size; i++){
        delete[] vec[i];
    }
    delete[] vec;
}

double Molecule::compute_ccsd(Matrix C, Matrix F, Matrix isqrt_S){
    
    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;
    Eigen::SelfAdjointEigenSolver<Matrix> solver(Fp);
    Matrix evals = solver.eigenvalues();

    int nso = 2 * norb;
    
    // convert from AO spatial to MO spatial
    double *moeri = spatial_atom(eri, C);
    // convert from MO spatial to MO spin
    double ****mospin = spatial_to_spin(moeri);
    
    Matrix Fs = Matrix::Zero(nso, nso);

    int noso = 2 * calc_nomo();

    for (int p = 0; p < nso; p++){
        for (int q = 0; q < nso; q++){
            if (p == q) {
                Fs(p, q) = evals(p/2);
            }
        }
    }

    double ****D_ijab = create4dmat(nso);
    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    D_ijab[i][j][a][b] = Fs(i, i) + Fs(j, j) - Fs(a, a) - Fs(b, b);
                }
            }
        }
    }

    double **D_ia = create2dvec(nso);
    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            D_ia[i][a] = Fs(i, i) - Fs(a, a);
        }
    }

    double ****t_ijab = create4dmat(nso);
    
    // initialize t_ijab
    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){                 
                    t_ijab[i][j][a][b] = mospin[i][j][a][b] / (Fs(i, i) + Fs(j, j) - Fs(a, a) - Fs(b, b));
                }
            }
        }
    }

    // initialize t_ia to zeros
    double **t_ia = create2dvec(nso);

    Matrix Fae = Matrix::Zero(nso, nso);
    Matrix Fmi = Matrix::Zero(nso, nso);
    Matrix Fme = Matrix::Zero(nso, nso);

    double ****Wmnij = create4dmat(nso);
    double ****Wabef = create4dmat(nso); 
    double ****Wmbej = create4dmat(nso);
    
    double E_curr = calc_ccsd_energy(Fs, mospin, t_ijab, t_ia);
    double E_prev = 0.0;
    double delta_E = E_prev - E_curr;

    int count = 0;
    printf("%20.12f\n", calc_ccsd_energy(Fs, mospin, t_ijab, t_ia));

    while (count < MAXITER && abs(delta_E) > DELTA_1){
        
        update_interm(mospin, t_ijab, t_ia, Fs, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej);
        updateT(mospin, t_ijab, t_ia, Fs, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej, D_ijab, D_ia);
        E_prev = E_curr;
        E_curr = calc_ccsd_energy(Fs, mospin, t_ijab, t_ia);
        delta_E = E_prev - E_curr;
        count++;
        printf("Iter = %2d Ecc = %20.12f\n", count, E_curr);
    }

    if (count == MAXITER){
        printf("CC iterations failed to converge\n");
        return 0.0;
    } else {
        printf("CC iterations converged\n");
        printf("Ecc = %21.12f\n", E_curr);
    }

    free4dmat(mospin, nso);
    
    free4dmat(t_ijab, nso);
    free2dvec(t_ia, nso);
    free4dmat(D_ijab, nso);
    free2dvec(D_ia, nso);

    free4dmat(Wmnij, nso);
    free4dmat(Wabef, nso);
    free4dmat(Wmbej, nso);

    return E_curr;
}

double Molecule::calc_ccsd_energy(Matrix Fs, double ****mospin, double ****t_ijab, double **t_ia){
    double energy = 0.0;
    int noso = 2 * calc_nomo();
    int nso = 2 * norb;
    
    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            energy += Fs(i, a) * t_ia[i][a];
        }
    }

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    energy += mospin[i][j][a][b] * t_ijab[i][j][a][b] / 4;
                    energy += mospin[i][j][a][b] * t_ia[i][a] * t_ia[j][b] / 2;
                }
            }
        }
    }
    return energy;
}

double Molecule::tau(double ****t_ijab, double **t_ia, int i, int j, int a, int b){
    return t_ijab[i][j][a][b] + t_ia[i][a] * t_ia[j][b] - t_ia[i][b] * t_ia[j][a];
}

double Molecule::taut(double ****t_ijab, double **t_ia, int i, int j, int a, int b){
    return t_ijab[i][j][a][b] + (t_ia[i][a] * t_ia[j][b] - t_ia[i][b] * t_ia[j][a]) / 2;
}

void Molecule::update_interm(double ****mospin, double ****t_ijab, double **t_ia, 
                             Matrix Fs, Matrix &Fae, Matrix &Fmi, Matrix &Fme,
                             double ****Wmnij, double ****Wabef, double ****Wmbej){
    
    int noso = 2 * calc_nomo();
    int nso = 2 * norb;
    
    // update Fae matrix (Equation 3)
    for (int a = noso; a < nso; a++){
        for (int e = noso; e < nso; e++){
            // Term 1
            Fae(a, e) = (1 - (a == e)) * Fs(a, e);
            for (int m = 0; m < noso; m++){
                // Term 2
                Fae(a, e) -= Fs(m, e) * t_ia[m][a] / 2;
                
                for (int f = noso; f < nso; f++){
                    // Term 3
                    Fae(a, e) += t_ia[m][f] * mospin[m][a][f][e];
                    
                    for (int n = 0; n < noso; n++){
                        // Term 4
                        Fae(a, e) -= taut(t_ijab, t_ia, m, n, a, f) * mospin[m][n][e][f] / 2;
                    }
                }
            }
        }
    }
    
    // update Fmi matrix (Equation 4)
    for (int m = 0; m < noso; m++){
        for (int i = 0; i < noso; i++){
            // Term 1
            Fmi(m, i) = (1 - (m == i)) * Fs(m, i);
            for (int e = noso; e < nso; e++){
                // Term 2
                Fmi(m, i) += t_ia[i][e] * Fs(m, e) / 2;
                for (int n = 0; n < noso; n++){
                    // Term 3
                    Fmi(m, i) += t_ia[n][e] * mospin[m][n][i][e];    
                    for (int f = noso; f < nso; f++){
                        // Term 4
                        Fmi(m, i) += taut(t_ijab, t_ia, i, n, e, f) * mospin[m][n][e][f] / 2;
                    }
                }
            }
        }
    }

    // update Fme matrix (Equation 5)
    for (int m = 0; m < noso; m++){
        for (int e = noso; e < nso; e++){
            // Term 1
            Fme(m, e) = Fs(m, e);
            for (int n = 0; n < noso; n++){
                for (int f = noso; f < nso; f++){
                    // Term 2
                    Fme(m, e) += t_ia[n][f] * mospin[m][n][e][f];
                }
            }
        }
    }

    // update Wmnij matrix (Equation 6)
    for (int m = 0; m < noso; m++){
        for (int n = 0; n < noso; n++){
            for (int i = 0; i < noso; i++){
                for (int j = 0; j < noso; j++){
                    // Term 1
                    Wmnij[m][n][i][j] = mospin[m][n][i][j];
                    for (int e = noso; e < nso; e++){
                        // Term 2
                        Wmnij[m][n][i][j] += t_ia[j][e] * mospin[m][n][i][e] - 
                                             t_ia[i][e] * mospin[m][n][j][e];
                        for (int f = noso; f < nso; f++){
                            // Term 3
                            Wmnij[m][n][i][j] += tau(t_ijab, t_ia, i, j, e, f) * 
                                                 mospin[m][n][e][f] / 4;
                        }
                    }
                }
            }
        }
    }

    // update Wabef matrix (Equation 7)
    for (int a = noso; a < nso; a++){
        for (int b = noso; b < nso; b++){
            for (int e = noso; e < nso; e++){
                for (int f = noso; f < nso; f++){
                    // Term 1
                    Wabef[a][b][e][f] = mospin[a][b][e][f];
                    for (int m = 0; m < noso; m++){
                        // Term 2
                        Wabef[a][b][e][f] -= t_ia[m][b] * mospin[a][m][e][f] -
                                             t_ia[m][a] * mospin[b][m][e][f];
                        for (int n = 0; n < noso; n++){
                            // Term 3
                            Wabef[a][b][e][f] += tau(t_ijab, t_ia, m, n, a, b) * 
                                                 mospin[m][n][e][f] / 4;
                        }
                    }
                }
            }
        }
    }

    // update Wmbej matrix (Equation 8)
    for (int m = 0; m < noso; m++){
        for (int b = noso; b < nso; b++){
            for (int e = noso; e < nso; e++){
                for (int j = 0; j < noso; j++){
                    // Term 1
                    Wmbej[m][b][e][j] = mospin[m][b][e][j];
                    for (int f = noso; f < nso; f++){
                        // Term 2
                        Wmbej[m][b][e][j] += t_ia[j][f] * mospin[m][b][e][f];
                    }
                    
                    for (int n = 0; n < noso; n++){
                        // Term 3
                        Wmbej[m][b][e][j] -= t_ia[n][b] * mospin[m][n][e][j];
                        
                        for (int f = noso; f < nso; f++){
                            // Term 4
                            Wmbej[m][b][e][j] -= (t_ijab[j][n][f][b] / 2 + t_ia[j][f] * t_ia[n][b]) * 
                                                  mospin[m][n][e][f];
                        }
                    }
                }
            }
        }
    }
}

void Molecule::updateT(double ****mospin, double ****t_ijab, double **t_ia, Matrix Fs, Matrix Fae, 
                       Matrix Fmi, Matrix Fme, double ****Wmnij, double ****Wabef, double ****Wmbej,
                       double ****D_ijab, double **D_ia){
    
    int noso = 2 * calc_nomo();
    int nso = 2 * norb;

    double **tmp2d = create2dvec(nso);

    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            // Term 1
            tmp2d[i][a] = Fs(i, a);
            for (int e = noso; e < nso; e++){
                // Term 2
                tmp2d[i][a] += t_ia[i][e] * Fae(a, e);
            }
            
            for (int m = 0; m < noso; m++){
                // Term 3
                tmp2d[i][a] += t_ia[m][a] * Fmi(m, i);
                for (int e = noso; e < nso; e++){
                    // Term 4
                    tmp2d[i][a] += t_ijab[i][m][a][e] * Fme(m, e);

                    for (int f = noso; f < nso; f++){
                        // Term 6
                        tmp2d[i][a] -= t_ijab[i][m][e][f] * mospin[m][a][e][f] / 2;
                    }

                    for (int n = 0; n < noso; n++){
                        // Term 7
                        tmp2d[i][a] -= t_ijab[m][n][a][e] * mospin[n][m][e][i] / 2;
                    }
                }
            }
    
            for (int n = 0; n < noso; n++){
                for (int f = noso; f < nso; f++){
                    // Term 5
                    tmp2d[i][a] -= t_ia[n][f] * mospin[n][a][i][f];
                }
            }

            t_ia[i][a] = tmp2d[i][a] / D_ia[i][a];
        }   
    }
    free2dvec(tmp2d, nso);

    double ****tmp4d = create4dmat(nso);

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    // Term 1
                    tmp4d[i][j][a][b] = mospin[i][j][a][b];
                    for (int e = noso; e < nso; e++){
                        // Term 2.1
                        tmp4d[i][j][a][b] += t_ijab[i][j][a][e] * Fae(b, e) - 
                                             t_ijab[i][j][b][e] * Fae(a, e);
                        // Term 7
                        tmp4d[i][j][a][b] += t_ia[i][e] * mospin[a][b][e][j] - 
                                             t_ia[j][e] * mospin[a][b][e][i];
                        for (int m = 0; m < noso; m++){
                            // Term 2.2
                            tmp4d[i][j][a][b] -= (t_ijab[i][j][a][e] * t_ia[m][b] * Fme(m, e) - 
                                                  t_ijab[i][j][b][e] * t_ia[m][a] * Fme(m, e)) / 2;
                        }

                        for (int f = noso; f < nso; f++){
                            // Term 5
                            tmp4d[i][j][a][b] += tau(t_ijab, t_ia, i, j, e, f) * Wabef[a][b][e][f] / 2;
                        }
                    }

                    for (int m = 0; m < noso; m++){
                        // Term 3.1
                        tmp4d[i][j][a][b] -= t_ijab[i][m][a][b] * Fmi(m, j) - 
                                             t_ijab[j][m][a][b] * Fmi(m, i);
                        
                        // Term 8
                        tmp4d[i][j][a][b] -= t_ia[m][a] * mospin[m][b][i][j] - 
                                             t_ia[m][b] * mospin[m][a][i][j];
                        
                        for (int e = noso; e < nso; e++){
                            // Term 3.2
                            tmp4d[i][j][a][b] -= (t_ijab[i][m][a][b] * t_ia[j][e] * Fme(m, e) -
                                                  t_ijab[j][m][a][b] * t_ia[i][e] * Fme(m, e)) / 2;

                            // Term 6 
                            tmp4d[i][j][a][b] += t_ijab[i][m][a][e] * Wmbej[m][b][e][j] - t_ia[i][e] * t_ia[m][a] * mospin[m][b][e][j];
                            tmp4d[i][j][a][b] -= t_ijab[j][m][a][e] * Wmbej[m][b][e][i] - t_ia[j][e] * t_ia[m][a] * mospin[m][b][e][i];
                            tmp4d[i][j][a][b] -= t_ijab[i][m][b][e] * Wmbej[m][a][e][j] - t_ia[i][e] * t_ia[m][b] * mospin[m][a][e][j];
                            tmp4d[i][j][a][b] += t_ijab[j][m][b][e] * Wmbej[m][a][e][i] - t_ia[j][e] * t_ia[m][b] * mospin[m][a][e][i];
                        }

                        for (int n = 0; n < noso; n++){
                            // Term 4
                            tmp4d[i][j][a][b] += tau(t_ijab, t_ia, m, n, a, b) * Wmnij[m][n][i][j] / 2;
                        }
                    }
                    t_ijab[i][j][a][b] = tmp4d[i][j][a][b] / D_ijab[i][j][a][b];
                }
                
            }
        }
    }

    free4dmat(tmp4d, nso);
}
