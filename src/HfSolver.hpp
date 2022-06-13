#ifndef HFSOLVER_H
#define HFSOLVER_H

#include "type.h"
#include "EnergySolver.hpp"

class HfSolver : public EnergySolver{
    std::string dir;
    int norb = 0, natom = 0;
    Atom *atoms;
    double enuc = 0.0;
    int *ioff;
    double **s, **t, **v, *eri, **ham;
    double **mux, **muy, **muz;
    Matrix C, D, F, isqrt_S;

    bool toprint;

    HfSolver(Molecule &m, bool toprint);
    ~HfSolver();
    void read_one_electron();
    void read_two_electron();
    double** readMatrix(std::string path);
    void hamiltonian();
    int calc_norb(std::string path);
    void lookupTable();
    void read_dipole();
    double calc_hf_energy();
    void initialize();
    int calc_nomo();
    void updateFock();
    void updateDensity(Matrix &new_D);
    double compute();
    void compute_dipole();
};

#endif /* HFSOLVER_H */