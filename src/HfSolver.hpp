#ifndef HFSOLVER_H
#define HFSOLVER_H

#include "type.h"
#include "EnergySolver.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"

class HfSolver : public EnergySolver{
public:
    
    HfSolver(Molecule &m, bool toprint);
    double compute();

protected:
    Matrix C, D, F, isqrt_S;

    double* spatial_atom();
    double**** spatial_to_spin(double *moeri);
    Matrix get_eval();

private:
    bool computed;

    double calc_hf_energy() const;
    void initialize();
    void updateFock();
    void updateDensity(Matrix &new_D);
    void compute_dipole() const;
};

#endif /* HFSOLVER_H */