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
    // Helper::print_matrix(e);
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

void DIIS::extrap(Matrix &ext) const{
    int size = (count > MAXERR) ? MAXERR : count;

    Matrix B(size + 1, size + 1);
    build_B(B);
    Eigen::VectorXd b(B.cols());
    for (int i = 0; i < B.cols(); i++){
        b[i] = 0;
    }
    b[B.cols() - 1] = -1;
    Eigen::VectorXd c = B.colPivHouseholderQr().solve(b);
    
    // std::cout << c << std::endl;

    for (int i = 0; i < c.size() - 1; i++){
        ext += c[i] * mats[i];
    }
}

void DIIS::build_B(Matrix &B) const{

    for (int i = 0; i < B.rows() - 1; i++){
        for (int j = 0; j <= i; j++){
            B(i, j) = (err[i].transpose() * err[j])(0, 0);
            B(j, i) = B(i, j);  //symmetry
        }
    }

    for (int i = 0; i < B.rows() - 1; i++){
        B(B.rows() - 1, i) = -1.0;
        B(i, B.cols() - 1) = -1.0;
    }

    B(B.rows() - 1, B.cols() - 1) = 0.0;
    // Helper::print_matrix(B);
}