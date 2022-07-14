#ifndef DIIS_H
#define DIIS_H

#include "type.h"

class DIIS{

public:
    int dim;
    
    DIIS();
    void add(const Matrix &mat, const Matrix &e);
    void extrap(Matrix &ext) const;
    double calc_rms(Matrix &e) const;
    
private:
    int count;
    std::vector<Matrix> err;
    std::vector<Matrix> mats;
    void build_B(Matrix &B) const;
    void shift();
};

#endif /* DIIS_H */
