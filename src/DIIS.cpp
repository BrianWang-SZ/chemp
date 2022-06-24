#include "DIIS.hpp"
#include "type.h"
#include "Helper.hpp"
#include <iostream>

#define MAXERR 8

DIIS::DIIS(){
    err = new Matrix[MAXERR];
    mats = new Matrix[MAXERR];
    count = 0;
}

void DIIS::add(Matrix &mat, Matrix &e){
    Helper::print_matrix(e);
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

void DIIS::extrap(Matrix &ext) const{
    Matrix B = build_B();
    Eigen::VectorXd b(B.cols());
    for (int i = 0; i < B.cols(); i++){
        b[i] = 0;
    }
    b[B.cols() - 1] = -1;
    Eigen::VectorXd c = B.householderQr().solve(b);
    std::cout << c << std::endl;
    Matrix ext(mats[0].rows(), mats[0].cols());

    for (int i = 0; i < c.size() - 1; i++){
        ext += c[i] * mats[i];
    }
}

Matrix DIIS::build_B() const{
    int size = (count > MAXERR) ? MAXERR : count;
    Matrix B(size + 1, size + 1);
    for (int i = 0; i < size; i++){
        for (int j = 0; j <= i; j++){
            B(i, j) = (err[i].transpose() * err[j]).trace();
            B(j, i) = B(i, j);  //symmetry
        }
    }

    for (int i = 0; i < size; i++){
        B(size, i) = -1;
        B(i, size) = -1;
    }
    B(size, size) = 0.0;

    return B;
}
