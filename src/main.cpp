#include "Molecule.hpp"
#include "GeomSolver.hpp"
#include "VibSolver.hpp"

int main (int argc, char **argv){
    const char *s = "../input"; 
    Molecule m(s);
    GeomSolver gs(m);
    gs.solve();
    VibSolver vs(m);
    vs.solve();
}
