#ifndef HFSOLVER_H
#define HFSOLVER_H

#include "type.h"
#include "EnergySolver.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"

class HfSolver : public EnergySolver{
public:
    
    HfSolver(Molecule &m, bool toprint);
    ~HfSolver();
    double compute();

protected:
    Matrix C, D, F, isqrt_S;

    double* spatial_atom();

private:
    bool toprint;
    bool computed;

    void read_dipole(double **mux, double **muy, double **muz);
    double calc_hf_energy();
    void initialize();
    void updateFock();
    void updateDensity(Matrix &new_D);
    void compute_dipole();
};

#endif /* HFSOLVER_H */