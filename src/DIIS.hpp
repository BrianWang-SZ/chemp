#ifndef DIIS_H
#define DIIS_H

#include "type.h"

class DIIS{

public:
    DIIS();
    void add(Matrix &mat, Matrix &e);
    Matrix extrap();

private:
    int count;
    std::vector<Matrix> err;
    std::vector<Matrix> mats;
    Matrix build_B();
    void shift();
    
};

#endif /* DIIS_H */