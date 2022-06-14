#ifndef HFSOLVER_H
#define HFSOLVER_H

#include "type.h"
#include "EnergySolver.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"

class HfSolver : public EnergySolver{
public:
    std::string dir;
    int norb = 0, natom = 0, nomo = 0;
    Atom *atoms;
    double enuc = 0.0;
    int *ioff;
    double *eri;
    double **s, **t, **v, **ham;
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
    void updateFock();
    void updateDensity(Matrix &new_D);
    double compute();
    void compute_dipole();
};

#endif /* HFSOLVER_H */