#include "HfSolver.hpp"
#include "Mp2Solver.hpp"

#define INDEX(i,j) ((i>j) ? (((ioff[i])+(j)) : ((ioff[j])+(i))))

Mp2Solver::Mp2Solver(Molecule &m):
    HfSolver(m, false){
    moeri = spatial_atom();
}

Mp2Solver::~Mp2Solver(){
    delete[] moeri;
}

double Mp2Solver::compute(){
    Matrix evals = get_eval();

    double *moeri = spatial_atom();

    double E_mp2 = 0.0;

    for (int i = 0; i < nomo; i++){
        for (int a = nomo; a < norb; a++){
            for (int j = 0; j < nomo; j++){
                for (int b = nomo; b < norb; b++){
                    int ia = INDEX(i, a);
                    int jb = INDEX(j, b);
                    int ib = INDEX(i, b);
                    int ja = INDEX(j, a);

                    E_mp2 += moeri[INDEX(ia, jb)] * (2 * moeri[INDEX(ia, jb)] - 
                                                        moeri[INDEX(ib, ja)]) /
                             (evals(i) + evals(j) - evals(a) - evals(b));
                }
            }
        }
    }

    double E_scf = HfSolver::compute();
    printf("Escf = %20.12f\n"
            "Emp2 = %20.12f\n"
            "Etot = %20.12f\n", E_scf + enuc, E_mp2, E_scf + enuc + E_mp2);
    
    return E_mp2;
}