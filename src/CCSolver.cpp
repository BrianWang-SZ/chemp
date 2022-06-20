#include "CCSolver.hpp"
#include "type.h"
#include "Helper.hpp"
#include "HfSolver.hpp"
#include "Molecule.hpp"

#define INDEX(i, j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))
#define MAXITER 2
#define DELTA_1 1e-12

CCSolver::CCSolver(Molecule &m): 
    HfSolver(m, false), Fs(nso, nso){

    // convert from AO spatial to MO spatial
    moeri = spatial_atom();
    // convert from MO spatial to MO spin
    mospin = spatial_to_spin(moeri);

    initialize_Fs();
    initialize_D();
    initialize_T();
}

CCSolver::~CCSolver(){
    delete[] moeri;
    Helper::free4d(mospin, nso);
    
    Helper::free2d(t_ai, nso);
    Helper::free4d(t_abij, nso);
    Helper::free2d(D_ai, nso);
    Helper::free4d(D_abij, nso);
}

double CCSolver::calc_ccsd_energy(){
    double energy = 0.0;
    
    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            energy += Fs(i, a) * t_ai[a][i];
        }
    }

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    energy += mospin[i][j][a][b] * t_abij[a][b][i][j] / 4;
                    energy += mospin[i][j][a][b] * t_ai[a][i] * t_ai[b][j] / 2;
                }
            }
        }
    }
    return energy;
}

double CCSolver::compute(){
    Matrix Fae = Matrix::Zero(nso, nso);
    Matrix Fmi = Matrix::Zero(nso, nso);
    Matrix Fme = Matrix::Zero(nso, nso);

    double ****Wmnij = Helper::create4d(nso);
    double ****Wabef = Helper::create4d(nso); 
    double ****Wmbej = Helper::create4d(nso);
    
    double E_curr = calc_ccsd_energy();
    printf("%20.12f\n", E_curr);
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

void CCSolver::initialize_Fs(){
    Matrix evals = get_eval();

    for (int p = 0; p < nso; p++){
        for (int q = 0; q < nso; q++){
            if (p == q) {
                Fs(p, q) = evals(p/2);
            }
        }
    }
}

void CCSolver::initialize_D(){

    D_abij = Helper::create4d(nso);

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    D_abij[a][b][i][j] = Fs(i, i) + Fs(j, j) - Fs(a, a) - Fs(b, b);
                }
            }
        }
    }

    D_ai = Helper::create2d(nso);
    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            D_ai[a][i] = Fs(i, i) - Fs(a, a);
        }
    }
}

void CCSolver::initialize_T(){

    t_abij = Helper::create4d(nso);
    
    // initialize t_ijab
    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){                 
                    t_abij[a][b][i][j] = mospin[i][j][a][b] / (Fs(i, i) + Fs(j, j) - Fs(a, a) - Fs(b, b));
                }
            }
        }
    }

    t_ai = Helper::create2d(nso);
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
                Fae(a, e) -= Fs(m, e) * t_ai[a][m] / 2;
                
                for (int f = noso; f < nso; f++){
                    // Term 3
                    Fae(a, e) += t_ai[f][m] * mospin[m][a][f][e];
                    
                    for (int n = 0; n < noso; n++){
                        // Term 4
                        Fae(a, e) -= taut(a, f, m, n) * mospin[m][n][e][f] / 2;
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
                Fmi(m, i) += t_ai[e][i] * Fs(m, e) / 2;
                for (int n = 0; n < noso; n++){
                    // Term 3
                    Fmi(m, i) += t_ai[e][n] * mospin[m][n][i][e];    
                    for (int f = noso; f < nso; f++){
                        // Term 4
                        Fmi(m, i) += taut(e, f, i, n) * mospin[m][n][e][f] / 2;
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
                    Fme(m, e) += t_ai[f][n] * mospin[m][n][e][f];
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
                        Wmnij[m][n][i][j] += t_ai[e][j] * mospin[m][n][i][e] - 
                                             t_ai[e][i] * mospin[m][n][j][e];
                        for (int f = noso; f < nso; f++){
                            // Term 3
                            Wmnij[m][n][i][j] += tau(e, f, i, j) * 
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
                        Wabef[a][b][e][f] -= t_ai[b][m] * mospin[a][m][e][f] -
                                             t_ai[a][m] * mospin[b][m][e][f];
                        for (int n = 0; n < noso; n++){
                            // Term 3
                            Wabef[a][b][e][f] += tau(a, b, m, n) * 
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
                        Wmbej[m][b][e][j] += t_ai[f][j] * mospin[m][b][e][f];
                    }
                    
                    for (int n = 0; n < noso; n++){
                        // Term 3
                        Wmbej[m][b][e][j] -= t_ai[b][n] * mospin[m][n][e][j];
                        
                        for (int f = noso; f < nso; f++){
                            // Term 4
                            Wmbej[m][b][e][j] -= (t_abij[f][b][j][n] / 2 + t_ai[f][j] * t_ai[b][n]) * 
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
            tmp2d[a][i] = Fs(i, a);
            for (int e = noso; e < nso; e++){
                // Term 2
                tmp2d[a][i] += t_ai[e][i] * Fae(a, e);
            }
            
            for (int m = 0; m < noso; m++){
                // Term 3
                tmp2d[a][i] += t_ai[a][m] * Fmi(m, i);
                for (int e = noso; e < nso; e++){
                    // Term 4
                    tmp2d[a][i] += t_abij[a][e][i][m] * Fme(m, e);

                    for (int f = noso; f < nso; f++){
                        // Term 6
                        tmp2d[a][i] -= t_abij[e][f][i][m] * mospin[m][a][e][f] / 2;
                    }

                    for (int n = 0; n < noso; n++){
                        // Term 7
                        tmp2d[a][i] -= t_abij[a][e][m][n] * mospin[n][m][e][i] / 2;
                    }
                }
            }
    
            for (int n = 0; n < noso; n++){
                for (int f = noso; f < nso; f++){
                    // Term 5
                    tmp2d[a][i] -= t_ai[f][n] * mospin[n][a][i][f];
                }
            }
        }   
    }
    

    // update T2 matrix (Equation 2)
    double ****tmp4d = Helper::create4d(nso);

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    // Term 1
                    tmp4d[a][b][i][j] = mospin[i][j][a][b];
                    for (int e = noso; e < nso; e++){
                        //Term 2.1
                        tmp4d[a][b][i][j] += t_abij[a][e][i][j] * Fae(b, e) - 
                                             t_abij[b][e][i][j] * Fae(a, e);
                        // Term 7
                        tmp4d[a][b][i][j] += t_ai[e][i] * mospin[a][b][e][j] - 
                                             t_ai[e][j] * mospin[a][b][e][i];
                        for (int m = 0; m < noso; m++){
                            // Term 2.2
                            // tmp4d[a][b][i][j] -= (t_abij[a][e][i][j] * t_ai[b][m] * Fme(m, e) - 
                            //                       t_abij[b][e][i][j] * t_ai[a][m] * Fme(m, e)) / 2;
                        }

                        for (int f = noso; f < nso; f++){
                            // Term 5
                            tmp4d[a][b][i][j] += tau(e, f, i, j) * Wabef[a][b][e][f] / 2;
                        }
                    }

                    for (int m = 0; m < noso; m++){
                        // Term 3.1
                        tmp4d[a][b][i][j] -= t_abij[a][b][i][m] * Fmi(m, j) - 
                                             t_abij[a][b][j][m] * Fmi(m, i);
                        
                        // Term 8
                        tmp4d[a][b][i][j] -= t_ai[a][m] * mospin[m][b][i][j] - 
                                             t_ai[b][m] * mospin[m][a][i][j];
                        
                        for (int e = noso; e < nso; e++){
                            // Term 3.2
                            tmp4d[a][b][i][j] -= (t_abij[a][b][i][m] * t_ai[e][j] * Fme(m, e) -
                                                  t_abij[a][b][j][m] * t_ai[e][i] * Fme(m, e)) / 2;

                            // Term 6 
                            tmp4d[a][b][i][j] += t_abij[a][e][i][m] * Wmbej[m][b][e][j] - t_ai[e][i] * t_ai[a][m] * mospin[m][b][e][j];
                            tmp4d[a][b][i][j] -= t_abij[a][e][j][m] * Wmbej[m][b][e][i] - t_ai[e][j] * t_ai[a][m] * mospin[m][b][e][i];
                            tmp4d[a][b][i][j] -= t_abij[b][e][i][m] * Wmbej[m][a][e][j] - t_ai[e][i] * t_ai[b][m] * mospin[m][a][e][j];
                            tmp4d[a][b][i][j] += t_abij[b][e][j][m] * Wmbej[m][a][e][i] - t_ai[e][j] * t_ai[b][m] * mospin[m][a][e][i];
                        }

                        for (int n = 0; n < noso; n++){
                            // Term 4
                            tmp4d[a][b][i][j] += tau(a, b, m, n) * Wmnij[m][n][i][j] / 2;
                        }
                    }
                } 
            }
        }
    }

    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            t_ai[a][i] = tmp2d[a][i] / D_ai[a][i];
        }
    }
    Helper::free2d(tmp2d, nso);

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    t_abij[a][b][i][j] = tmp4d[a][b][i][j] / D_abij[a][b][i][j];
                }
            }
        }
    }
    Helper::free4d(tmp4d, nso);
}

double CCSolver::tau(int a, int b, int i, int j){
    return t_abij[a][b][i][j] + t_ai[a][i] * t_ai[b][j] - t_ai[b][i] * t_ai[a][j];
}

double CCSolver::taut(int a, int b, int i, int j){
    return t_abij[a][b][i][j] + (t_ai[a][i] * t_ai[b][j] - t_ai[b][i] * t_ai[a][j]) / 2;
}