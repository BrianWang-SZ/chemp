#include "Molecule.hpp"
#include "GeomSolver.hpp"

int main (int argc, char **argv){
    const char *s = "../input/h2o/DZ"; 
    Molecule m(s);
    GeomSolver gs(m);
    gs.solve();
}
