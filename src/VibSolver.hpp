#ifndef VIBSOLVER_H
#define VIBSOLVER_H

#include "Atom.hpp"
#include "type.h"
#include "Molecule.hpp"

class VibSolver{
public:
    int natom = 0;
    Atom *atoms;
    int size;
    double **hes;

    VibSolver(Molecule &m);
    ~VibSolver();

    void read_hes(std::string dir);
    void weight();
    void compute();
    void solve();
};
#endif /* VIBSOLVER_H */