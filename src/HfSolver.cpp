#include "HfSolver.hpp"
#include "Molecule.hpp"
#include "Helper.hpp"
#include "EnergySolver.hpp"
#include "DIIS.hpp"

#define MAXITER 100
#define DELTA_1 1e-12
#define DELTA_2 1e-11

 #define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

HfSolver::HfSolver(Molecule &m, bool toprint):
    EnergySolver(m, false), C(norb, norb), D(norb, norb),
    F(norb, norb), isqrt_S(norb, norb){
    this -> toprint = toprint;
    computed = false;
    if (!toprint) compute();
}

Matrix HfSolver::get_eval(){
    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;
    Eigen::SelfAdjointEigenSolver<Matrix> solver(Fp);
    return solver.eigenvalues();
}

void HfSolver::compute_dipole(){
    double **mux, **muy, **muz;
    
    std::string path = dir + "/mux.dat";
    mux = readMatrix(path.c_str());

    path = dir + "/muy.dat";
    muy = readMatrix(path.c_str());

    path = dir + "/muz.dat";
    muz = readMatrix(path.c_str());

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

    Helper::free2d(mux, norb);
    Helper::free2d(muy, norb);
    Helper::free2d(muz, norb);
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
    if (computed && !toprint) return calc_hf_energy();
    
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

    DIIS d;

    Matrix S(norb, norb);

    for (int i = 0; i < norb; i++){
        for (int j = 0; j < norb; j++){
            S(i, j) = s[i][j];
        }
    }

    while (count < MAXITER && (abs(delta_E) >= DELTA_1 || rms >= DELTA_2)){

        E_prev = E_curr;
        
        Matrix new_D(norb, norb);

        /* DIIS optimization starts*/
        if(count >= 2){
            Matrix oldF = F;
            F = d.extrap();
            updateFock();

            Matrix e = F * D * S - S * D * F;
            Helper::print_matrix(e);

            updateDensity(new_D);
            F = oldF;
            
            d.add(F, e);
            
        } else {
            updateFock();
            
            if(toprint && count == 0){
                printf("\tFock Matrix:\n\n");
                Helper::print_matrix(F);
            }

            updateDensity(new_D);
            
            Matrix e = F * D * S - S * D * F;
            Helper::print_matrix(e);
            d.add(F, e);
            
        }
        /* DIIS optimization ends*/

        /****/
        // updateFock();

        // Matrix e = F * D * S - S * D * F;
        // Helper::print_matrix(e);
            
        // if(toprint && count == 0){
        //     printf("\tFock Matrix:\n\n");
        //     Helper::print_matrix(F);
        // }
        
        // updateDensity(new_D);
        /****/

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

    // initialize S matrix
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

    isqrt_S = evecs_S * isqrt_A * evecs_S.transpose();

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

double* HfSolver::spatial_atom(){
    
    int max = INDEX(norb, norb);
    double *moeri = new double[INDEX(max, max) + 1];

    for (int i = 0; i < INDEX(max, max); i++){
        moeri[i] = 0.0;
    }

    double ****M = Helper::create4d(norb);
    double ****N = Helper::create4d(norb);
    double ****P = Helper::create4d(norb);
    double ****Q = Helper::create4d(norb);
    
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

    Helper::free4d(M, norb);
    Helper::free4d(N, norb);
    Helper::free4d(P, norb);
    Helper::free4d(Q, norb);

    return moeri;
}

double**** HfSolver::spatial_to_spin(double *moeri){

    double ****mospin = Helper::create4d(nso);

    for (int p = 0; p < nso; p++){
        for (int q = 0; q < nso; q++){
            for (int r = 0; r < nso; r++){
                for (int s = 0; s < nso; s++){
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
