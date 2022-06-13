
#include "Atom.hpp"
#include "type.h"

class VibSolver{
public:
    int natom = 0;
    Atom *atoms;
    int size;
    Vec2d hes;

    VibSolver(const char *, const char*);
    void VibSolver::read_hes(const char* path);
    void VibSolver::weight();
    void VibSolver::compute();
    void VibSolver::solve();
};