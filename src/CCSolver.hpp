#ifndef CCSOLVER_H
#define CCSOLVER_H

#include "EnergySolver.hpp"
#include "Molecule.hpp"
#include "HfSolver.hpp" 

class CCSolver : public HfSolver{
public:
    int nso, noso;
    double *moeri;
    double ****mospin;
    double **t_ia, **D_ia;
    double ****t_ijab, ****D_ijab;
    Matrix Fs;

    CCSolver(Molecule &m);
    ~CCSolver();
    double compute();
    void spatial_atom();
    void spatial_to_spin();
    void initialize_Fs();
    void initialize_D();
    void initialize_T();

    double tau(int i, int j, int a, int b);
    double taut(int i, int j, int a, int b);
    
    double calc_ccsd_energy();
    
    void update_interm(Matrix &Fae, Matrix &Fmi, Matrix &Fme,
                       double ****Wmnij, double ****Wabef, double ****Wmbej);

    void updateT( Matrix Fae, Matrix Fmi, Matrix Fme, 
                  double ****Wmnij, double ****Wabef, double ****Wmbej);

};

#endif /* CCSOLVER_H */