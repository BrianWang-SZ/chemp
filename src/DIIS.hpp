#ifndef DIIS_H
#define DIIS_H

#include "type.h"

class DIIS{

public:
    DIIS();
    void add(Matrix &mat, Matrix &e);
    void extrap(Matrix &ext) const;

private:
    int count;
    Matrix *err;
    Matrix *mats;
    Matrix build_B() const;
    void shift();
    
};

#endif /* DIIS_H */