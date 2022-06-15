#ifndef MP2SOLVER_H
#define MP2SOLVER_H

class Mp2Solver : public HfSolver{
public:

    Mp2Solver(Molecule &m);
    ~Mp2Solver();
    double compute();

private:
    double *moeri;
}

#endif /* MP2SOLVER_H