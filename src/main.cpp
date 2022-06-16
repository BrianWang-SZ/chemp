#include "Molecule.hpp"
#include "GeomSolver.hpp"
#include "VibSolver.hpp"
#include "HfSolver.hpp"
#include "CCSolver.hpp"
#include "Mp2Solver.hpp"

int main (int argc, char **argv){
    const char *s = "../input/h2o/STO-3G"; 
    Molecule m(s);
    HfSolver hfs(m, true);
    hfs.compute();
    // GeomSolver gs(m);
    // gs.solve();
    // VibSolver vs(m);
    // vs.solve();
    //Mp2Solver ms(m);
    //ms.compute();
    //CCSolver ccs(m);
    //ccs.compute();
    
}
