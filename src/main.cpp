#include "Molecule.hpp"
#include "GeomSolver.hpp"
#include "VibSolver.hpp"

int main (int argc, char **argv){
    const char *s = "../input/h2o_geom.txt"; 
    Molecule m(s);
    GeomSolver gs(m);
    gs.solve();
    VibSolver vs(m, "../input/h2o_hessian.txt");
    vs.solve();
}
