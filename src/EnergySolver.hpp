#ifndef ENERGYSOLVER_H
#define ENERGYSOLVER_H

#include <string>
#include "type.h"
#include "Molecule.hpp"

class EnergySolver{

public:
    EnergySolver(Molecule &m, bool toprint);
    ~EnergySolver();
    virtual double compute() = 0;

protected:
    
    bool toprint;
    
    std::string dir;
    int natom;

    int norb, nomo,
        nso, noso;

    Atom *atoms;
    double enuc;
    int *ioff;
    double *eri;
    double **s, **t, **v, **ham;

    double** readMatrix(std::string path) const;

private:
    void read_one_electron();
    void read_two_electron();
    void hamiltonian();
    int calc_norb(std::string path) const;
    void lookupTable();
};

#endif /* ENERGYSOLVER_H */