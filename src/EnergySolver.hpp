#ifndef ENERGYSOLVER_H
#define ENERGYSOLVER_H

#include <string>
#include "type.h"
#include "Molecule.hpp"

class EnergySolver{

public:
    EnergySolver(Molecule &m);
    ~EnergySolver();
    virtual double compute() = 0;

protected:
    std::string dir;
    int natom = 0;

    int norb = 0, nomo = 0,
        nso = 0, noso = 0;
        
    Atom *atoms;
    double enuc = 0.0;
    int *ioff;
    double *eri;
    double **s, **t, **v, **ham;

    double* spatial_atom();
    double**** spatial_to_spin(double *moeri);
    double** readMatrix(std::string path);

private:
    void read_one_electron();
    void read_two_electron();
    void hamiltonian();
    int calc_norb(std::string path);
    void lookupTable();
};

#endif /* ENERGYSOLVER_H */