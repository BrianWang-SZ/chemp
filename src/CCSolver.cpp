#include "CCSolver.hpp"
#include "type.h"
#include "Helper.hpp"
#include "HFSolver.hpp"
#include "Molecule.hpp"

#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))
#define MAXITER 100
#define DELTA_1 1e-12

CCSolver::CCSolver(Molecule &m): 
    HFSolver(m,toprint=false){

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
    Helper::free2d(Fs, nso);
    
    Helper::free2d(t_ia, nso);
    Helper::free4d(t_ijab, nso);
    Helper::free2d(D_ia, nso);
    Helper::free4d(D_ijab, nso);
}

double CCSolver::calc_ccsd_energy() const{
    double energy = 0.0;
    
    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            energy += Fs[i][a] * t_ia[i][a];
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

double CCSolver::compute(){
    double **Fae = Helper::create2d(nso);
    double **Fmi = Helper::create2d(nso);
    double **Fme = Helper::create2d(nso);

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

    Helper::free2d(Fae, nso);
    Helper::free2d(Fmi, nso);
    Helper::free2d(Fme, nso);

    Helper::free4d(Wmnij, nso);
    Helper::free4d(Wabef, nso);
    Helper::free4d(Wmbej, nso);

    return E_curr;
}

void CCSolver::initialize_Fs(){
    Matrix evals = get_eval();

    Fs = Helper::create2d(nso);

    for (int p = 0; p < nso; p++){
        for (int q = 0; q < nso; q++){
            if (p == q) {
                Fs[p][q] = evals(p / 2);
            }
        }
    }
}

void CCSolver::initialize_D(){

    // Equation 13
    D_ijab = Helper::create4d(nso);
    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    D_ijab[i][j][a][b] = Fs[i][i] + Fs[j][j] - Fs[a][a] - Fs[b][b];
                }
            }
        }
    }

    // Equation 12
    D_ia = Helper::create2d(nso);
    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            D_ia[i][a] = Fs[i][i] - Fs[a][a];
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
                    t_ijab[i][j][a][b] = mospin[i][j][a][b] / (Fs[i][i] + Fs[j][j] - Fs[a][a] - Fs[b][b]);
                }
            }
        }
    }

    t_ia = Helper::create2d(nso);
}

void CCSolver::update_interm(double **Fae, double **Fmi, double **Fme,
                             double ****Wmnij, double ****Wabef, double ****Wmbej) const{
    
    // update Fae matrix (Equation 3)
    for (int a = noso; a < nso; a++){
        for (int e = noso; e < nso; e++){
            // Term 1
            Fae[a][e] = (1 - (a == e)) * Fs[a][e];
            
            // Term 2
            for (int m = 0; m < noso; m++){
                
                Fae[a][e] -= Fs[m][e] * t_ia[m][a] / 2;
            }

            // Term 3
            for (int m = 0; m < noso; m++){    
                for (int f = noso; f < nso; f++){
                    Fae[a][e] += t_ia[m][f] * mospin[m][a][f][e];
                }
            }

            // Term 4
            for (int m = 0; m < noso; m++){ 
                for (int n = 0; n < noso; n++){
                    for (int f = noso; f < nso; f++){   
                        Fae[a][e] -= taut(m, n, a, f) * mospin[m][n][e][f] / 2;
                    }
                }
            }
        }
    }
    
    // update Fmi matrix (Equation 4)
    for (int m = 0; m < noso; m++){
        for (int i = 0; i < noso; i++){
            // Term 1
            Fmi[m][i] = (1 - (m == i)) * Fs[m][i];
            
            // Term 2
            for (int e = noso; e < nso; e++){
                Fmi[m][i] += t_ia[i][e] * Fs[m][e] / 2;
            }

            // Term 3
            for (int e = noso; e < nso; e++){
                for (int n = 0; n < noso; n++){
                    Fmi[m][i] += t_ia[n][e] * mospin[m][n][i][e];    
                }
            }

            // Term 4
            for (int n = 0; n < noso; n++){
                for (int e = noso; e < nso; e++){
                    for (int f = noso; f < nso; f++){
                        Fmi[m][i] += taut(i, n, e, f) * mospin[m][n][e][f] / 2;
                    }
                }
            }
        }
    }

    // update Fme matrix (Equation 5)
    for (int m = 0; m < noso; m++){
        for (int e = noso; e < nso; e++){
            // Term 1
            Fme[m][e] = Fs[m][e];

            // Term 2
            for (int n = 0; n < noso; n++){
                for (int f = noso; f < nso; f++){
                    Fme[m][e] += t_ia[n][f] * mospin[m][n][e][f];
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

                    // Term 2
                    for (int e = noso; e < nso; e++){
                        Wmnij[m][n][i][j] += t_ia[j][e] * mospin[m][n][i][e] - 
                                             t_ia[i][e] * mospin[m][n][j][e];
                    }

                    // Term 3
                    for (int e = noso; e < nso; e++){    
                        for (int f = noso; f < nso; f++){
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

                    // Term 2
                    for (int m = 0; m < noso; m++){
                        Wabef[a][b][e][f] -= t_ia[m][b] * mospin[a][m][e][f] -
                                             t_ia[m][a] * mospin[b][m][e][f];
                    }

                    // Term 3
                    for (int m = 0; m < noso; m++){
                        for (int n = 0; n < noso; n++){
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

                    // Term 2
                    for (int f = noso; f < nso; f++){
                        Wmbej[m][b][e][j] += t_ia[j][f] * mospin[m][b][e][f];
                    }

                    // Term 3
                    for (int n = 0; n < noso; n++){
                        Wmbej[m][b][e][j] -= t_ia[n][b] * mospin[m][n][e][j];
                    }

                    // Term 4
                    for (int n = 0; n < noso; n++){
                        for (int f = noso; f < nso; f++){
                            Wmbej[m][b][e][j] -= (t_ijab[j][n][f][b] / 2 + t_ia[j][f] * t_ia[n][b]) * 
                                                  mospin[m][n][e][f];
                        }
                    }
                }
            }
        }
    }
}

void CCSolver::updateT(double **Fae, double **Fmi, double **Fme, 
                       double ****Wmnij, double ****Wabef, double ****Wmbej){

    // update T1 matrix (Equation 1)
    double **tmp2d = Helper::create2d(nso);

    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            // Term 1
            tmp2d[i][a] = Fs[i][a];

            // Term 2
            for (int e = noso; e < nso; e++){
                tmp2d[i][a] += t_ia[i][e] * Fae[a][e];
            }

            // Term 3
            for (int m = 0; m < noso; m++){
                tmp2d[i][a] += t_ia[m][a] * Fmi[m][i];
            }

            // Term 4
            for (int m = 0; m < noso; m++){
                for (int e = noso; e < nso; e++){
                    
                    tmp2d[i][a] += t_ijab[i][m][a][e] * Fme[m][e];
                }
            }

            // Term 5
            for (int n = 0; n < noso; n++){
                for (int f = noso; f < nso; f++){
                    tmp2d[i][a] -= t_ia[n][f] * mospin[n][a][i][f];
                }
            }

            // Term 6
            for (int m = 0; m < noso; m++){
                for (int e = noso; e < nso; e++){
                    for (int f = noso; f < nso; f++){
                        tmp2d[i][a] -= t_ijab[i][m][e][f] * mospin[m][a][e][f] / 2;
                    }
                }
            }

            // Term 7
            for (int m = 0; m < noso; m++){
                for (int e = noso; e < nso; e++){
                    for (int n = 0; n < noso; n++){
                        tmp2d[i][a] -= t_ijab[m][n][a][e] * mospin[n][m][e][i] / 2;
                    }
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
                    tmp4d[i][j][a][b] = mospin[i][j][a][b];
                    
                    // Term 2.1
                    for (int e = noso; e < nso; e++){
                        tmp4d[i][j][a][b] += t_ijab[i][j][a][e] * Fae[b][e] - 
                                             t_ijab[i][j][b][e] * Fae[a][e];
                    }

                    // Term 2.2
                    for (int e = noso; e < nso; e++){
                        for (int m = 0; m < noso; m++){
                            tmp4d[i][j][a][b] -= (t_ijab[i][j][a][e] * t_ia[m][b] * Fme[m][e] - 
                                                  t_ijab[i][j][b][e] * t_ia[m][a] * Fme[m][e]) / 2;
                        }
                    }

                    // Term 3.1
                    for (int m = 0; m < noso; m++){
                        tmp4d[i][j][a][b] -= t_ijab[i][m][a][b] * Fmi[m][j] - 
                                             t_ijab[j][m][a][b] * Fmi[m][i];
                    }
                    
                    // Term 3.2
                    for (int m = 0; m < noso; m++){
                        for (int e = noso; e < nso; e++){
                            tmp4d[i][j][a][b] -= (t_ijab[i][m][a][b] * t_ia[j][e] * Fme[m][e] -
                                                  t_ijab[j][m][a][b] * t_ia[i][e] * Fme[m][e]) / 2;
                        }
                    }
                    
                    // Term 4
                    for (int m = 0; m < noso; m++){
                        for (int n = 0; n < noso; n++){
                            tmp4d[i][j][a][b] += tau(m, n, a, b) * Wmnij[m][n][i][j] / 2;
                        }
                    }

                    // Term 5
                    for (int e = noso; e < nso; e++){
                        for (int f = noso; f < nso; f++){
                            tmp4d[i][j][a][b] += tau(i, j, e, f) * Wabef[a][b][e][f] / 2;
                        }
                    }

                    // Term 6 
                    for (int m = 0; m < noso; m++){
                        for (int e = noso; e < nso; e++){
                            tmp4d[i][j][a][b] += t_ijab[i][m][a][e] * Wmbej[m][b][e][j] -
                                                 t_ia[i][e] * t_ia[m][a] * mospin[m][b][e][j];
                            
                            tmp4d[i][j][a][b] -= t_ijab[j][m][a][e] * Wmbej[m][b][e][i] -
                                                 t_ia[j][e] * t_ia[m][a] * mospin[m][b][e][i];
                            
                            tmp4d[i][j][a][b] -= t_ijab[i][m][b][e] * Wmbej[m][a][e][j] -
                                                 t_ia[i][e] * t_ia[m][b] * mospin[m][a][e][j];
                            
                            tmp4d[i][j][a][b] += t_ijab[j][m][b][e] * Wmbej[m][a][e][i] -
                                                 t_ia[j][e] * t_ia[m][b] * mospin[m][a][e][i];
                        }
                    }
                    
                    // Term 7
                    for (int e = noso; e < nso; e++){
                        tmp4d[i][j][a][b] += t_ia[i][e] * mospin[a][b][e][j] - 
                                             t_ia[j][e] * mospin[a][b][e][i];
                    }
                    
                    // Term 8
                    for (int m = 0; m < noso; m++){
                        tmp4d[i][j][a][b] -= t_ia[m][a] * mospin[m][b][i][j] - 
                                             t_ia[m][b] * mospin[m][a][i][j];
                    }
                } 
            }
        }
    }

    for (int i = 0; i < noso; i++){
        for (int a = noso; a < nso; a++){
            if (D_ia[i][a] == 0){
                t_ia[i][a] = 0;
	    } else {
		t_ia[i][a] = tmp2d[i][a] / D_ia[i][a];
            }
	}
    }
    Helper::free2d(tmp2d, nso);

    for (int i = 0; i < noso; i++){
        for (int j = 0; j < noso; j++){
            for (int a = noso; a < nso; a++){
                for (int b = noso; b < nso; b++){
                    if (D_ijab[i][j][a][b] == 0){
	                t_ijab[i][j][a][b] = 0;
		    } else {
		        t_ijab[i][j][a][b] = tmp4d[i][j][a][b] / D_ijab[i][j][a][b];
                    }
		}
            }
        }
    }
    Helper::free4d(tmp4d, nso);
}

// Equation 10
double CCSolver::tau(int i, int j, int a, int b) const{
    return t_ijab[i][j][a][b] + t_ia[i][a] * t_ia[j][b] - t_ia[i][b] * t_ia[j][a];
}

// Equation 9
double CCSolver::taut(int i, int j, int a, int b) const{
    return t_ijab[i][j][a][b] + (t_ia[i][a] * t_ia[j][b] - t_ia[i][b] * t_ia[j][a]) / 2;
}
