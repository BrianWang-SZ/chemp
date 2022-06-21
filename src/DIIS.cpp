#include "DIIS.hpp"
#include "type.h"
#include "Helper.hpp"

#define MAXERR 8

DIIS::DIIS():
    err(MAXERR), mats(MAXERR){
    count = 0;
}

void DIIS::add(Matrix mat, Matrix e){
    

    printf("here\n");
    e.resize(e.rows() * e.cols(), 1);
    
    if (count >= MAXERR) {
        shift();
        err[MAXERR - 1] = e;
        mats[MAXERR - 1] = mat;
    } else {
        err[count] = e;
        mats[count] = mat;
    }

    count++;
}

void DIIS::shift(){
    for (int i = 0; i < MAXERR - 1; i++){
        err[i] = err[i + 1];
        mats[i] = mats[i + 1];
    }
}

Matrix DIIS::extrap(){
    Matrix B = build_B();
    Eigen::VectorXd b(B.cols());
    b[B.cols() - 1] = -1;
    Eigen::VectorXd c = B.householderQr().solve(b);
    Matrix ext(mats[0].rows(), mats[0].cols());

    for (int i = 0; i < c.size() - 1; i++){
        ext += c[i] * mats[i];
    }
    return ext;
}

Matrix DIIS::build_B(){
    int size = (count > MAXERR) ? MAXERR : count;
    Matrix B(size + 1, size + 1);
    for (int i = 0; i < size; i++){
        for (int j = 0; j <= i; j++){
            B(i, j) = (err[i].transpose() * err[j])(0, 0);
            B(j, i) = B(i, j);  //symmetry
        }
    }

    for (int i = 0; i < size; i++){
        B(size, i) = -1;
        B(i, size) = -1;
    }

    return B;
}
