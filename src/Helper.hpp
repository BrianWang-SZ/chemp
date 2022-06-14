#ifndef HELPER_H
#define HELPER_H

#include "type.h"

namespace Helper {
    void print_matrix(double **mat, int size);
    void print_matrix(Matrix mat);
    double** create2d(int size);
    double calc_rms(Matrix mat, Matrix newmat);
    void free2d(double **vec, int size);
    double**** create4d(int size);
    void free4d(double ****mat, int size);
}

#endif /* HELPER_H */