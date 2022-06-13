#include "Atom.hpp"
#include "Molecule.hpp"
#include "GeomSolver.hpp"

#include "const.h"
#include "masses.h"
#include "math.h"

// float comparison precision
#define EPSILON 1e-4

GeomSolver::GeomSolver(Molecule &molec){
    natom = molec.natom;
    atoms = molec.atoms;
    bond_length();
    unitvec();
    bond_angle();
}

GeomSolver::~GeomSolver(){

    // free length
    for (int i = 0; i < natom; i++){
        delete[] length[i];
    }
    delete[] length;

    // free angle
    for (int i = 0; i < natom; i++){
        for (int j = 0; j < natom; j++){
            delete[] angle[i][j];
        }
        delete[] angle[i];
    }

    delete[] angle;
    // free vec
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < natom; j++){
            delete[] vec[i][j];
        }
        delete[] vec[i];
    }

    delete[] vec;
}

void GeomSolver::print_coord(){
    printf("Number of atoms: %d\n", natom);

    printf("Input Cartesian Coordinates:\n");
    for (int i = 0; i < natom; i++){
        printf("%d%19.12f%19.12f%19.12f\n", atoms[i].zval, atoms[i].x, atoms[i].y, atoms[i].z);
    }
}

void GeomSolver::bond_length(){
    
    // allocate memory for 2D array
    double **len = new double*[natom];
    for (int i = 0; i < natom; i++){
        len[i] = new double[natom];
    }

    for (int i = 0; i < natom; i++){
        for (int j = 0; j < natom; j++){
            Atom a = atoms[i];
            Atom b = atoms[j];
            len[i][j] = sqrt(pow(a.x - b.x, 2) +
                             pow(a.y - b.y, 2) +
                             pow(a.z - b.z, 2));
        }
    }

    length = len;
}

void GeomSolver::print_len(){
    printf("\nInteratomic distances (bohr):\n");
    for (int i = 1; i < natom; i++){
        for (int j = 0; j < i; j++){
            printf("%2d%2d%9.5f\n", i, j, length[i][j]);
        }
    }
    printf("\n");
}

void GeomSolver::unitvec(){
    vec = new double**[3];

    for (int i = 0; i < 3; i++){
        vec[i] = new double*[natom];
        for (int j = 0; j < natom; j++){
            vec[i][j] = new double[natom];
        }
    }

    // calculate unit vector
    for (int i = 0; i < natom; i++){
        for (int j = 0; j < natom; j++){
            
            if (i != j){
                double leng = length[i][j];
            
                vec[0][i][j] = -(atoms[i].x - atoms[j].x) / leng;
                vec[1][i][j] = -(atoms[i].y - atoms[j].y) / leng;
                vec[2][i][j] = -(atoms[i].z - atoms[j].z) / leng;
            }
    
        }
    }
}

void GeomSolver::print_vec(){
    for (int i = 0; i < natom; i++){
        for (int j = 0; j < natom; j++){
            printf("%f ", vec[0][i][j]);    
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < natom; i++){
        for (int j = 0; j < natom; j++){
            printf("%f ", vec[1][i][j]);    
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < natom; i++){
        for (int j = 0; j < natom; j++){
            printf("%f ", vec[2][i][j]);    
        }
        printf("\n");
    }
}

void GeomSolver::bond_angle(){

    // 3d array initialization
    double ***angle = new double**[natom];
    for (int i = 0; i < natom; i++){
        angle[i] = new double*[natom];

        for (int j = 0; j < natom; j++){
            angle[i][j] = new double[natom];
        }
    }

    for (int i = 0; i < natom; i++){
        for (int j = 0; j < natom; j++){
            for (int k = 0; k < natom; k++){
                // eliminate duplicate atoms
                if (i != j && j != k && i != k){
                    angle[i][j][k] = acos(vec[0][j][i] * vec[0][j][k] +
                                          vec[1][j][i] * vec[1][j][k] +
                                          vec[2][j][i] * vec[2][j][k]); 
                }
            }
        }
    }

    this -> angle = angle;
}

void GeomSolver::print_angle(){
    printf("Bond angles:\n");
    for (int i = 0; i < natom; i++){
        for (int j = i + 1; j < natom; j++){
            for (int k = j + 1; k < natom; k++){
                // eliminate duplicate atoms
                if (i != j && j != k && i != k){

                    // only consider proximate atoms
                    if (length[i][j] < 4.0 && length[j][k] < 4.0){
                        printf("%2d-%2d-%2d %11.6f\n", i, j, k, angle[i][j][k] / M_PI * 180);
                    } else if (length[i][j] < 4.0 && length[i][k] < 4.0){
                        printf("%2d-%2d-%2d %11.6f\n", j, i, k, angle[j][i][k] / M_PI * 180);
                    } else if (length[i][k] < 4.0 && length[j][k] < 4.0){
                        printf("%2d-%2d-%2d %11.6f\n", i, k, j, angle[i][k][j] / M_PI * 180);
                    }
                }
            }
        }
    }
    printf("\n"); 
}

void GeomSolver::oop_angle(){
    printf("Out-of-plane angles:\n\n");

    // for all possible central atoms
    for (int k = 0; k < natom; k++){

        // for all possible combinations of the other three atoms
        for (int i = 0; i < natom; i++){
            for (int j = i + 1; j < natom; j++){
                for (int l = j + 1; l < natom; l++){
                    
                    // eliminate duplicate atoms and distant atoms
                    if (i != j && i != k && i != l && j != k && 
                        j != l && k != l && length[i][k] < 4.0 && 
                        length[k][j] < 4.0 && length[k][l] < 4.0){
                        
                        // consider three possible oop atoms
                        oop_helper(i, j, k, l);
                        oop_helper(j, i, k, l);
                        oop_helper(l, i, k, j);                                            
                    }  
                }
            }
        }
    }
    printf("\n");
}

void GeomSolver::oop_helper(int i, int j, int k, int l){
    double prodx = vec[1][k][j] * vec[2][k][l] - 
                   vec[2][k][j] * vec[1][k][l];
    double prody = vec[2][k][j] * vec[0][k][l] - 
                   vec[0][k][j] * vec[2][k][l];
    double prodz = vec[0][k][j] * vec[1][k][l] - 
                   vec[1][k][j] * vec[0][k][l];

    double x = prodx * vec[0][k][i];
    double y = prody * vec[1][k][i];
    double z = prodz * vec[2][k][i];
    double result = (x + y + z) / sin(angle[j][k][l]);
        
    if (result < -1) result = -1;
    else if (result > 1) result = 1;

    printf("%2d-%2d-%2d-%2d%11.6f\n", i, j, k, l, asin(result) / M_PI * 180);   
}

void GeomSolver::tor_angle(){
    printf("Torsional angles:\n\n");
    for (int i = 0; i < natom; i++){
        for (int j = 0; j < i; j++){
            for (int k = 0; k < j; k++){
                for (int l = 0; l < k; l++){
                    if (i != j && i != k && i != l && j != k && 
                        j != l && k != l && length[i][j] < 4.0 && 
                        length[j][k] < 4.0 && length[k][l] < 4.0){
                        
                        double prodxl = vec[1][i][j] * vec[2][j][k] - 
                                        vec[2][i][j] * vec[1][j][k];
                        double prodyl = vec[2][i][j] * vec[0][j][k] - 
                                        vec[0][i][j] * vec[2][j][k];
                        double prodzl = vec[0][i][j] * vec[1][j][k] - 
                                        vec[1][i][j] * vec[0][j][k];

                        double prodxr = vec[1][j][k] * vec[2][k][l] - 
                                        vec[2][j][k] * vec[1][k][l];
                        double prodyr = vec[2][j][k] * vec[0][k][l] - 
                                        vec[0][j][k] * vec[2][k][l];
                        double prodzr = vec[0][j][k] * vec[1][k][l] - 
                                        vec[1][j][k] * vec[0][k][l];
                        
                        double result = (prodxl * prodxr + prodyl * prodyr + prodzl * prodzr) / 
                                         sin(angle[i][j][k]) / sin(angle[j][k][l]);
                        
                        if (result < -1) result = -1;
                        else if (result > 1) result = 1;
                        
                        printf("%2d-%2d-%2d-%2d%11.6f\n", i, j, k, l, acos(result) / M_PI * 180);
                    }
                }
            }
        }
    }
    printf("\n");
}

void GeomSolver::massctr(double *ctr){

    double tot = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    for (int i = 0; i < natom; i++){
        tot += an2masses[atoms[i].zval];
        x += an2masses[atoms[i].zval] * atoms[i].x;
        y += an2masses[atoms[i].zval] * atoms[i].y;
        z += an2masses[atoms[i].zval] * atoms[i].z;    
    }

    x /= tot;
    y /= tot;
    z /= tot;
    
    if (ctr != NULL){
        ctr[0] = x;
        ctr[1] = y;
        ctr[2] = z;

        return;
    }

    printf("Molecular center of mass:%11.8f%11.8f%11.8f\n", x, y, z);
}

void GeomSolver::moi(double *result){
   
    double ctr[3];
    massctr(ctr);

    // initialize tensor matrix
    Matrix tensor(3,3);

    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;

    for (int i = 0; i < natom; i++){

        // calculate relative position to mass center
        double x = atoms[i].x - ctr[0];
        double y = atoms[i].y - ctr[1];
        double z = atoms[i].z - ctr[2];

        double m = an2masses[atoms[i].zval];

        xx += m * (pow(y, 2) + pow(z, 2));
        yy += m * (pow(x, 2) + pow(z, 2));
        zz += m * (pow(x, 2) + pow(y, 2));

        xy += m * x * y;
        xz += m * x * z;
        yz += m * y * z;
    }

    tensor(0, 0) = xx;
    tensor(1, 1) = yy;
    tensor(2, 2) = zz;

    tensor(0, 1) = -xy;
    tensor(0, 2) = -xz;
    tensor(1, 2) = -yz;

    tensor(1, 0) = tensor(0, 1);
    tensor(2, 0) = tensor(0, 2);
    tensor(2, 1) = tensor(1, 2);

    Eigen::SelfAdjointEigenSolver<Matrix> solver(tensor);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();

    if (result == NULL){
        print_moi(tensor, evals);
    } else {
        for (int i = 0; i < 3; i++){
            result[i] = evals(i);
        }
    }
}

void GeomSolver::print_moi(Matrix tensor, Matrix evals){

    printf("\nMoment of inertia tensor:\n");

    printf("%12d%12d%12d\n", 1, 2, 3);

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            if (j == 0) printf("%4d", i + 1);
            printf("%12.7f", tensor(i, j));
        }
        printf("\n");
    }

    printf("\nPrincipal moments of inertia (amu * bohr^2):\n    ");
    for (double x : evals.reshaped()){
        printf("%12.6f", x);
    }

    printf("\n\nPrincipal moments of inertia (amu * aa^2):\n    ");
    for (double x : evals.reshaped()){
        printf("%12.6f", x * pow(BOHR / ANGSTROM, 2));
    }
  
    printf("\n\nPrincipal moments of inertia (g * cm^2):\n    ");
    for (double x : evals.reshaped()){
        printf("%13.6e", x * AMU * 1000 * pow(BOHR * 100, 2));
    }

    printf("\n\n");
    
    if (natom == 2) 
        printf("Molecule is diatomic.");
    else if (fabs(evals(0) < EPSILON))
        printf("Molecule is linear.");
    else if (fabs(evals(1) - evals(0)) < EPSILON && fabs(evals(2) - evals(1)) < EPSILON) 
        printf("Molecule is a spherical top.");
    else if (fabs(evals(1) - evals(0)) > EPSILON && fabs(evals(2) - evals(1)) < EPSILON)
        printf("Molecule is a prolate symmetric top.");
    else if (fabs(evals(1) - evals(0)) < EPSILON && fabs(evals(2) - evals(1)) > EPSILON)
        printf("Molecule is an oblate symmetric top.");
    else
        printf("Molecule is an asymmetric top.");
    
    printf("\n");
}

void GeomSolver::rotation(){
    double inert[3];
    moi(inert);
    double conv = (pow(BOHR, 2) * AMU);

    printf("\nRotational constants (MHz):\n");

    double a = PLANCK / (8 * pow(M_PI, 2) * LIGHT * inert[0] * conv) * LIGHT / 1e6;
    double b = PLANCK / (8 * pow(M_PI, 2) * LIGHT * inert[1] * conv) * LIGHT / 1e6;
    double c = PLANCK / (8 * pow(M_PI, 2) * LIGHT * inert[2] * conv) * LIGHT / 1e6;

    printf("\tA =%10.3f\t B =%10.3f\t C =%10.3f\n\n", a, b, c);

    printf("Rotational constants (cm-1):\n");

    a = PLANCK / (8 * pow(M_PI, 2) * LIGHT * inert[0] * conv) / 100;
    b = PLANCK / (8 * pow(M_PI, 2) * LIGHT * inert[1] * conv) / 100;
    c = PLANCK / (8 * pow(M_PI, 2) * LIGHT * inert[2] * conv) / 100;

    printf("\tA = %6.4f\t B = %6.4f\t C = %6.4f\n", a, b, c);
}

void GeomSolver::solve(){
    print_coord();
    print_len();
    print_angle();
    oop_angle();
    tor_angle();
    massctr(NULL);
    moi(NULL);
    rotation();
}