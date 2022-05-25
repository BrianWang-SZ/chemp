#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>
#include "masses.h"
#include "const.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include "molecule.hpp"

#define MAXLINE 200

int main(int argc, char **argv){
    const char *geom = "./input/3c1b_geom.txt";
    const char *hess = "./input/3c1b_hessian.txt";
    Molecule molec(geom, hess);
}

Molecule::Molecule(const char *geom, const char *hess){
    coord(geom);
    double **hes = hessian(hess);
    //print_matrix(hes);
    weight(hes);
    //print_matrix(hes);
    compute(hes);
    cleanup(hes);
}

void Molecule::coord(const char *path){

    // open input file
    FILE *in;
    
    if ((in = fopen(path, "r")) == NULL) {
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

double** Molecule::hessian(const char* path){
    // open input file
    FILE *in;
    if ((in = fopen(path, "r")) == NULL) {
        perror("fopen() error");
        exit(-1);
    }
    
    // read in the number of atoms
    int num;
    if (fscanf(in, " %d", &num) != 1) {
        
        perror("Error reading number of atoms\n");
        exit(-1);
    }

    if (num != natom){
        printf("Inconsistent geometry and matrix files!\n");
    }
    
    // number of row/col of hessian matrix
    this -> size = 3 * natom;

    // initialize hessian matrix
    double **hes = new double*[size];
    
    for (int i = 0; i < size; i++){
        hes[i] = new double[size];
    }

    int count = 0;

    // three entries per input line
    while (fscanf(in, "%lf %lf %lf", &hes[count / size][count % size], 
                                     &hes[count / size][count % size + 1],
                                     &hes[count / size][count % size + 2]) == 3){
        count +=3;
    }

    if (fclose(in)){
        perror("Error closing file stream!");
        exit(-1);
    }

    return hes;
}

void Molecule::weight(double **hes){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i / 3;
            int b = j / 3;
            hes[i][j] /= sqrt(an2masses[atoms[a].zval] * an2masses[atoms[b].zval]);
        }
    }
}

void Molecule::print_matrix(double **mat){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            printf("%6.3f", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void Molecule::compute(double **hes){
    Matrix H(size, size);

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            H(i, j) = hes[i][j];
        }
    }

    Eigen::SelfAdjointEigenSolver<Matrix> solver(H);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();

    printf("\nHessian eigenvalues (hartree/amu-bohr^2):\n");
    for (int i = size - 1; i >= 0; i--){
        printf("%4d%21.10f\n", i, evals(i));
    }

    double conv = HARTREE / (AMU * pow(BOHR, 2));
    printf("\nHarmonic vibrational frequencies (cm^-1):\n");
    for (int i = size - 1; i >= 0; i--){
        if (evals(i) < 0){
            std::complex<double> c(evals(i), 0);
            printf("%4d%10.4fi\n", i, (sqrt(c * conv) / (2 * M_PI * LIGHT) / 100.0).imag());
        } else {
            printf("%4d%10.4f\n", i, sqrt(evals(i) * conv) / (2 * M_PI * LIGHT) / 100);
        }
    }
    
}

void Molecule::cleanup(double **hes){
    
    for (int i = 0; i < size; i++){
        delete[] hes[i];
    }

    delete[] hes;
}
