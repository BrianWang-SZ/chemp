#ifndef HFSOLVER_H
#define HFSOLVER_H

#include "type.h"
#include "EnergySolver.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"

class HFSolver : public EnergySolver{
public:
    
    HFSolver(Molecule &m, bool toprint=true, bool useDIIS=true);
    double compute();

protected:
    Matrix C, D, F, isqrt_S;

    double* spatial_atom();
    double**** spatial_to_spin(double *moeri);
    Matrix get_eval();

private:
    bool useDIIS;
    double calc_hf_energy(Matrix &D, Matrix &F) const;
    void initialize();
    void updateFock(Matrix &F, const Matrix &D);
    void updateDensity(Matrix &new_D, const Matrix &F);
    void compute_dipole() const;
};

#endif /* HFSOLVER_H */
