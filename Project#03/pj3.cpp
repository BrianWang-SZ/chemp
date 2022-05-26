#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include "molecule.hpp"

#define MAXORB 100
#define MAXITER 64
#define DELTA_1 1e-12
#define DELTA_2 1e-11

#define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i)

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
    double **ham = new double*[norb];
    
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            ham[i] = new double[norb];
        }
    }

    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            ham[i][j] = t[i][j] + v[i][j];
        }
    }
    this -> ham = ham;

    printf("\tCore Hamiltonian:\n\n");
    print_matrix(ham);

}

void Molecule::read_one_electron(const char *dir){

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
    print_matrix(s);
    
    path = direc + "/t.dat";
    this -> t = readMatrix(path.c_str());
    printf("\tKinetic-Energy Integrals:\n\n");
    print_matrix(t);

    path = direc + "/v.dat";
    this -> v = readMatrix(path.c_str());
    printf("\tNuclear Attraction Integrals:\n\n");
    print_matrix(v);
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

    double val;

    int ij = INDEX(i, j);
    int kl = INDEX(k, l);
    int ijkl = INDEX(ij, kl);

    double *eri = new double[ijkl];

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
    double **mat = new double*[norb];
    
    for (int i = 0; i < norb; i++){
        mat[i] = new double[norb];
    }

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
    coord(dir);
    read_one_electron(dir);
    hamiltonian();
    read_two_electron(dir);
    read_dipole(dir);
    compute_HF();
}

int main (int argc, char **argv){
    setbuf(stdout, NULL);
    const char *s = "./input/h2o/STO-3G";
    Molecule m = Molecule(s);
}

void Molecule::print_matrix(double **mat){
    for (int i = 0; i < norb; i++){
        printf("%12d", i + 1);
    }

    printf("\n\n");

    for (int i = 0; i < norb; i++){
        printf("%5d", i + 1);
        for (int j = 0; j < norb; j++){
            printf("%12.7f", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void Molecule::print_matrix(Matrix mat){
    for (int i = 0; i < norb; i++){
        printf("%12d", i + 1);
    }

    printf("\n\n");

    for (int i = 0; i < norb; i++){
        printf("%5d", i + 1);
        for (int j = 0; j < norb; j++){
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

void Molecule::compute_HF(){

    Matrix S(norb, norb);

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

    Matrix isqrt_S = evecs_S * isqrt_A * evecs_S.transpose();

    printf("\tS^-1/2 Matrix:\n\n");
    print_matrix(isqrt_S);

    Matrix H(norb, norb);

    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            H(i, j) = ham[i][j];
        }
    }

    Matrix F = isqrt_S.transpose() * H * isqrt_S;

    printf("\tInitial F' Matrix:\n\n");
    print_matrix(F);

    solver.compute(F);
    Matrix evecs_F = solver.eigenvectors();
    Matrix evals_F = solver.eigenvalues();

    Matrix C = isqrt_S * evecs_F;

    printf("\tInitial C Matrix:\n\n");
    print_matrix(C);
    
    Matrix D(norb, norb);

    // calculate number of occupied orbitals
    int nelec = 0;
    for (int i = 0; i < natom; i++){
        nelec += atoms[i].zval;
    }
    int nomo = nelec / 2;

    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            double sum = 0.0;
            for (int k = 0; k < nomo; k++){
                sum += C(i, k) * C(j, k);
            }
            D(i, j) = sum;
        }
    }

    printf("\tInitial Density Matrix:\n\n");
    print_matrix(D);

    double E_elec_prev = 0.0;
    double E_elec_curr = calc_energy(D, H, F);

    double delta_E = E_elec_curr - E_elec_prev;
    double rms = 1.0;

    int count = 0;

    printf("Iter\t\tE(elec)  \t\tE(tot)  \t\tDelta(E)  \t\tRMS(D)\n");
    printf("%02d%21.12f%21.12f\n\n", count, E_elec_curr, E_elec_curr + enuc);

    while (delta_E > DELTA_1 || rms > DELTA_2 && count < MAXITER){
        // compute new Fock Matrix
        updateFock(F, H, D);

        if(count == 0){
            printf("\tFock Matrix:\n\n");
            print_matrix(F);
        }
        //compute new density matrix
        Matrix new_D(norb, norb);
        updateDensity(new_D, isqrt_S, F, C);

        rms = calc_rms(D, new_D);

        D = new_D;

        E_elec_prev = E_elec_curr;
        E_elec_curr = calc_energy(new_D, H, F);
        delta_E = E_elec_curr - E_elec_prev;

        count++;
        printf("%02d%21.12f%21.12f%21.12f%21.12f\n", count, E_elec_curr, 
                                                     E_elec_curr + enuc, delta_E, rms); 
    }
    
    print_matrix(D);
    //compute_dipole(D);
}

double Molecule::calc_rms(Matrix D, Matrix new_D){
    double rms;
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            rms += sqrt(pow(D(i, j) - new_D(i, j), 2));
        }
    }
    return rms;
}

double Molecule::calc_energy(Matrix D, Matrix H, Matrix F){
    double E = 0.0;
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            E += D(i, j) * (H(i, j) + F(i, j));
        }
    }
    return E;
}

void Molecule::updateFock(Matrix &F, Matrix H, Matrix D){
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            
            double sum = 0.0;
            
            for (int k = 0; k < norb; k++){
                for (int l = 0; l < norb; l++){
                    int ij = INDEX(i, j);
                    int kl = INDEX(k, l);
                    int ik = INDEX(i, k);
                    int jl = INDEX(j, l);

                    sum += D(k, l) * (2 * eri[INDEX(ij, kl)] - eri[INDEX(ik, jl)]);
                }
            }
            F(i, j) = H(i, j) + sum;
        }
    }
    cleanNan(F);
}

void Molecule::mobasis(Matrix C, Matrix F){
    Matrix moF(norb, norb);
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

void Molecule::updateDensity(Matrix &new_D, Matrix isqrt_S, Matrix F, Matrix C){
    
    F = isqrt_S.transpose() * F * isqrt_S;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(F);
    Matrix evecs_F = solver.eigenvectors();
    Matrix evals_F = solver.eigenvalues();

    cleanNan(evecs_F);

    C = isqrt_S * evecs_F;

    int nelec = 0;

    for (int i = 0; i < natom; i++){
        nelec += atoms[i].zval;
    }

    int nomo = nelec / 2;

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

void Molecule::cleanNan(Matrix &M){
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){            
            M(i, j) = (std::isnan(M(i, j))) ?  0 : M(i, j);
        }
    }
}

void Molecule::coord(const char *dir){

    std::string direc = dir;
    std::string path = direc + "/geom.dat";
    
    // open input file
    FILE *in;
    
    if ((in = fopen(path.c_str(), "r")) == NULL) {
        perror("fopen() error");
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

// FIXME
void Molecule::compute_dipole(Matrix D){
    double dx = 0.0, dy = 0.0, dz = 0.0;
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            dx += 2 * D(i, j) * mux[i][j];
            dy += 2 * D(i, j) * muy[i][j];
            dz += 2 * D(i, j) * muz[i][j];
        }
    }
    printf("Mux = %20.12f\nMuy = %20.12f\nMuz = %20.12f\n", dx, dy, dz);
}
