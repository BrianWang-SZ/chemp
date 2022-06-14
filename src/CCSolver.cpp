#include "CCSolver.hpp"
#include "type.h"
#include "Helper.hpp"
#include "HfSolver.hpp"
#include "Molecule.hpp"

#define INDEX(i, j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))
#define MAXITER 100
#define DELTA_1 1e-12

CCSolver::CCSolver(Molecule &m) : HfSolver(m, false){
    noso = 2 * nomo;
    nso = 2 * norb;

    HfSolver::compute();

    // convert from AO spatial to MO spatial
    spatial_atom();
    // convert from MO spatial to MO spin
    spatial_to_spin();
    
    initialize_Fs();
    initialize_D();
    initialize_T();
}

CCSolver::~CCSolver(){
    delete[] moeri;
    
    Helper::free4d(mospin, nso);
    
    Helper::free2d(t_ia, nso);
    Helper::free4d(t_ijab, nso);
    Helper::free2d(D_ia, nso);
    Helper::free4d(D_ijab, nso);
}

double CCSolver::compute(){
    Matrix Fae = Matrix::Zero(nso, nso);
    Matrix Fmi = Matrix::Zero(nso, nso);
    Matrix Fme = Matrix::Zero(nso, nso);

    double ****Wmnij = Helper::create4d(nso);
    double ****Wabef = Helper::create4d(nso); 
    double ****Wmbej = Helper::create4d(nso);
    
    double E_curr = calc_ccsd_energy();
    double E_prev = 0.0;
    double delta_E = E_prev - E_curr;

    int count = 0;

    while (count < MAXITER && abs(delta_E) > DELTA_1){
        
        update_interm(Fae, Fmi, Fme, Wmnij, Wabef, Wmbej);
        updateT(Fae, Fmi, Fme, Wmnij, Wabef, Wmbej);
        E_prev = E_curr;
        E_curr = calc_ccsd_energy();
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

    Helper::free4d(Wmnij, nso);
    Helper::free4d(Wabef, nso);
    Helper::free4d(Wmbej, nso);

    return E_curr;
}

void CCSolver::spatial_atom(){
    int max = INDEX(norb - 1, norb - 1);
    moeri = new double[INDEX(max, max) + 1];

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
}

void CCSolver::spatial_to_spin(){

    mospin = Helper::create4d(nso);

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
}

void CCSolver::initialize_Fs(){
    Matrix Fp = isqrt_S.transpose() * F * isqrt_S;
    Eigen::SelfAdjointEigenSolver<Matrix> solver(Fp);
    Matrix evals = solver.eigenvalues();

    Fs = Matrix::Zero(nso, nso);

    for (int p = 0; p < nso; p++){
        for (int q = 0; q < nso; q++){
            if (p == q) {
                Fs(p, q) = evals(p/2);
            }
        }
    }
}

void CCSolver::initialize_D(){

    D_ijab = Helper::create4d(nso);

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    D_ijab[i][j][a][b] = Fs(i, i) + Fs(j, j) - Fs(a, a) - Fs(b, b);
                }
            }
        }
    }

    D_ia = Helper::create2d(nso);
    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            D_ia[i][a] = Fs(i, i) - Fs(a, a);
        }
    }
}

void CCSolver::initialize_T(){

    t_ijab = Helper::create4d(nso);
    
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
    t_ia = Helper::create2d(nso);
}

void CCSolver::update_interm(Matrix &Fae, Matrix &Fmi, Matrix &Fme,
                             double ****Wmnij, double ****Wabef, double ****Wmbej){
    
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
                        Fae(a, e) -= taut(m, n, a, f) * mospin[m][n][e][f] / 2;
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
                        Fmi(m, i) += taut(i, n, e, f) * mospin[m][n][e][f] / 2;
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
                            Wmnij[m][n][i][j] += tau(i, j, e, f) * 
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
                            Wabef[a][b][e][f] += tau(m, n, a, b) * 
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

void CCSolver::updateT(Matrix Fae, Matrix Fmi, Matrix Fme, 
                       double ****Wmnij, double ****Wabef, double ****Wmbej){

    // update T1 matrix (Equation 1)
    double **tmp2d = Helper::create2d(nso);

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
    Helper::free2d(tmp2d, nso);

    // update T2 matrix (Equation 2)
    double ****tmp4d = Helper::create4d(nso);

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
                            tmp4d[i][j][a][b] += tau(i, j, e, f) * Wabef[a][b][e][f] / 2;
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
                            tmp4d[i][j][a][b] += tau(m, n, a, b) * Wmnij[m][n][i][j] / 2;
                        }
                    }
                    t_ijab[i][j][a][b] = tmp4d[i][j][a][b] / D_ijab[i][j][a][b];
                } 
            }
        }
    }
    Helper::free4d(tmp4d, nso);
}