#include "HfSolver.hpp"
#include "Molecule.hpp"
#include "Helper.hpp"
#include <fstream>
#include <algorithm>

#define MAXORB 100
#define MAXITER 100
#define DELTA_1 1e-12
#define DELTA_2 1e-11

#define INDEX(i, j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

HfSolver::HfSolver(Molecule &m, bool toprint){
    dir = m.dir;
    natom = m.natom;
    atoms = m.atoms;
    nomo = m.nomo;
    this -> toprint = toprint;
    
    read_one_electron();
    hamiltonian();
    lookupTable();
    read_two_electron();
    read_dipole();
}

HfSolver::~HfSolver(){
    delete[] eri;
    delete[] ioff;

    Helper::free2d(s, norb);
    Helper::free2d(v, norb);
    Helper::free2d(t, norb);
    
    Helper::free2d(ham, norb);
    Helper::free2d(mux, norb);
    Helper::free2d(muy, norb);
    Helper::free2d(muz, norb);
}

void HfSolver::read_two_electron(){
    FILE *in;
    std::string path = dir + "/eri.dat";
    
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

    double *eri = new double[ijkl + 1];

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

int HfSolver::calc_norb(std::string path){
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

double** HfSolver::readMatrix(std::string path){

    if (this -> norb == 0) {
        int norb = calc_norb(path);
        if (norb == -1){
            perror("Exceeds maximum orbitals!");
            exit(-1);
        }
        this -> norb = norb;
    }
    
    FILE *in;
    double **mat = Helper::create2d(norb);

    if ((in = fopen(path.c_str(), "r")) == NULL) {
        fprintf(stderr, "Error opening %s\n", path.c_str());
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
        fprintf(stderr, "Inconsistent input file %s\n", path.c_str());
        perror("");
        exit(-1);
    }

    if (fclose(in)){
        fprintf(stderr, "Error closing %s\n", path.c_str());
        perror("");
        exit(-1);
    }
    return mat;
}

void HfSolver::read_one_electron(){

    FILE *in;
    std::string path = dir + "/enuc.dat";
    
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

    path = dir + "/s.dat";
    this -> s = readMatrix(path);
    
    path = dir + "/t.dat";
    this -> t = readMatrix(path);

    path = dir + "/v.dat";
    this -> v = readMatrix(path);
    

    if (toprint){
        printf("\n\tOverlap Integrals:\n\n");
        Helper::print_matrix(s, norb);

        printf("\tKinetic-Energy Integrals:\n\n");
        Helper::print_matrix(t, norb);

        printf("\tNuclear Attraction Integrals:\n\n");
        Helper::print_matrix(v, norb);
    }
}

void HfSolver::hamiltonian(){
    double **ham = Helper::create2d(norb);
    
    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            ham[i][j] = t[i][j] + v[i][j];
        } 
    }
    this -> ham = ham;

    if (toprint){
        printf("\tCore Hamiltonian:\n\n");
        Helper::print_matrix(ham, norb);
    }
}

// lookup array
void HfSolver::lookupTable(){
    int* ioff = new int[MAXORB];

    ioff[0] = 0;

    for (int i = 1; i < MAXORB; i++){
        ioff[i] = ioff[i - 1] + i;
    }

    this -> ioff = ioff;
}

void HfSolver::read_dipole(){
    std::string path = dir + "/mux.dat";
    this -> mux = readMatrix(path.c_str());

    path = dir + "/muy.dat";
    this -> muy = readMatrix(path.c_str());

    path = dir + "/muz.dat";
    this -> muz = readMatrix(path.c_str());
}

void HfSolver::compute_dipole(){
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
    initialize();

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