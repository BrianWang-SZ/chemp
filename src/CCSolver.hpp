#ifndef CCSOLVER_H
#define CCSOLVER_H

#include "Molecule.hpp"
#include "HfSolver.hpp" 

class CCSolver : public HfSolver{
public:

    CCSolver(Molecule &m);
    ~CCSolver();
    double compute();

private:
    double *moeri;
    double ****mospin;
    
    double **t_ia, **D_ia;
    double ****t_ijab, ****D_ijab;
    double **Fs;
    
    void initialize_Fs();
    void initialize_D();
    void initialize_T();

    double tau(int i, int j, int a, int b);
    double taut(int i, int j, int a, int b);
    
    double calc_ccsd_energy();
    
    void update_interm(double **Fae, double **Fmi, double **Fme,
                       double ****Wmnij, double ****Wabef, double ****Wmbej);

    void updateT( double **Fae, double **Fmi, double **Fme, 
                  double ****Wmnij, double ****Wabef, double ****Wmbej);

};

#endif /* CCSOLVER_H */
