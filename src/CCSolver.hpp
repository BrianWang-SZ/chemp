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
    
    double **t_ai, **D_ai;
    double ****t_abij, ****D_abij;
    Matrix Fs;
    
    void initialize_Fs();
    void initialize_D();
    void initialize_T();

    double tau(int a, int b, int i, int j);
    double taut(int a, int b, int i, int j);
    
    double calc_ccsd_energy();
    
    void update_interm(Matrix &Fae, Matrix &Fmi, Matrix &Fme,
                       double ****Wmnij, double ****Wabef, double ****Wmbej);

    void updateT( Matrix Fae, Matrix Fmi, Matrix Fme, 
                  double ****Wmnij, double ****Wabef, double ****Wmbej);

};

#endif /* CCSOLVER_H */
