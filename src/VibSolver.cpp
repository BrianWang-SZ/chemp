#include "VibSolver.hpp"
#include "masses.h"
#include "Molecule.hpp"
#include "Atom.hpp"
#include "const.h"
#include <math.h>
#include "type.h"

VibSolver::VibSolver(Molecule &m, const char *hess){
    natom = m.natom;
    atoms = m.atoms;

    read_hes(hess);
}

VibSolver::~VibSolver(){
    for (int i = 0; i < size; i++){
	delete[] hes[i];
    }
    delete[] hes;
}

void VibSolver::read_hes(const char* path){
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
    hes = new double*[size];

    for (int i = 0; i < size; i++){
	hes[i] = new double[size];
    }

    int count = 0;

    // three entries per input line
    while (fscanf(in, "%lf %lf %lf", &hes[count / size][count % size], 
                                     &hes[count / size][count % size + 1],
                                     &hes[count / size][count % size + 2]) == 3){
        count += 3;
    }

    if (fclose(in)){
        perror("Error closing file stream!");
        exit(-1);
    }
}

void VibSolver::weight(){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i / 3;
            int b = j / 3;
            hes[i][j] /= sqrt(an2masses[atoms[a].zval] * 
			 an2masses[atoms[b].zval]);
        }
    }
}

void VibSolver::compute(){
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

void VibSolver::solve(){
    weight();
    compute();
}
