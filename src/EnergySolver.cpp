#include "EnergySolver.hpp"
#include "Helper.hpp"
#include "HfSolver.hpp"
#include <fstream>
#include <algorithm>

#define MAXORB 100
#define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i)

EnergySolver::EnergySolver(Molecule &m, bool toprint){
    this -> toprint = toprint;

    dir = m.dir;
    natom = m.natom;
    atoms = m.atoms;
    
    nomo = m.nomo;
    noso = 2 * nomo;
    
    read_one_electron();
    nso = 2 * norb;
    hamiltonian();
    
    lookupTable();
    read_two_electron();
}

EnergySolver::~EnergySolver(){
    delete[] eri;
    delete[] ioff;

    Helper::free2d(s, norb);
    Helper::free2d(v, norb);
    Helper::free2d(t, norb);
    Helper::free2d(ham, norb);
}

void EnergySolver::read_two_electron(){
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

int EnergySolver::calc_norb(std::string path){
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

double** EnergySolver::readMatrix(std::string path){

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

void EnergySolver::read_one_electron(){

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

    path = dir + "/s.dat";
    this -> s = readMatrix(path);
    
    path = dir + "/t.dat";
    this -> t = readMatrix(path);

    path = dir + "/v.dat";
    this -> v = readMatrix(path);
    

    if (toprint){
        printf("Nuclear repulsion energy = %20.10f\n", enuc);
        
        printf("\n\tOverlap Integrals:\n\n");
        Helper::print_matrix(s, norb);

        printf("\tKinetic-Energy Integrals:\n\n");
        Helper::print_matrix(t, norb);

        printf("\tNuclear Attraction Integrals:\n\n");
        Helper::print_matrix(v, norb);
    }
}

void EnergySolver::hamiltonian(){
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
void EnergySolver::lookupTable(){
    int* ioff = new int[MAXORB];

    ioff[0] = 0;

    for (int i = 1; i < MAXORB; i++){
        ioff[i] = ioff[i - 1] + i;
    }

    this -> ioff = ioff;
}