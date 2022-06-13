#include "Atom.hpp"
#include "type.h"
#include "Molecule.hpp"

class VibSolver{
public:
    int natom = 0;
    Atom *atoms;
    int size;
    double **hes;

    VibSolver(Molecule &m, const char*);
    ~VibSolver();

    void read_hes(const char* path);
    void weight();
    void compute();
    void solve();
};
