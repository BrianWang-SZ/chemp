#ifndef GEOMSOLVER_H
#define GEOMSOLVER_H

#include "Atom.hpp"
#include "Molecule.hpp"
#include "type.h"

class GeomSolver{
public:
    int natom;
    Atom *atoms;
    double **length;
    double ***angle;
    double ***vec;
    
    GeomSolver(Molecule &molec);
    ~GeomSolver();

    void bond_length();
    void unitvec();
    void bond_angle();
    void oop_angle();
    void tor_angle();
    void massctr(double*);
    void moi(double*);
    void rotation();
    
    void print_coord();
    void print_len();
    void print_angle();
    void print_vec();
    void print_moi(Matrix, Matrix);

    void solve();

private:
    void oop_helper(int i, int j, int k, int l);
};

#endif /* GEOMSOLVER_H */