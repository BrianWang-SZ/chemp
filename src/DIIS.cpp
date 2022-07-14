#include "DIIS.hpp"
#include "type.h"
#include "Helper.hpp"
#include <iostream>

#define MAXERR 6

DIIS::DIIS(){
    err = std::vector<Matrix>(MAXERR);
    mats = std::vector<Matrix>(MAXERR);
    count = 0;
}

void DIIS::add(const Matrix &mat, const Matrix &e){
    
    if (count >= MAXERR) shift();
    
    err[count] = e;
    mats[count] = mat;
    
    count++;
}

void DIIS::shift(){
    
    for (int i = 0; i < MAXERR - 1; i++){
        err[i] = err[i + 1];
        mats[i] = mats[i + 1];
    }
    count--;
}

void DIIS::extrap(Matrix &ext) const{
    if (count < 2){
        return;
    } else {
        ext.setZero();
        
        Matrix B = Matrix(count + 1, count + 1);
        build_B(B);
        
        Eigen::VectorXd b(count + 1);
        for (int i = 0; i < count; i++){
            b[i] = 0;
        }
        b[count] = -1;
        
        Eigen::VectorXd c = B.colPivHouseholderQr().solve(b);
        
        for (int i = 0; i < count; i++){
            ext += c[i] * mats[i];
        }
    }
}

void DIIS::build_B(Matrix &B) const{
    
    B.setZero();
    
    for (int i = 0; i < B.rows() - 1; i++){
        for (int j = 0; j <= i; j++){
            B(i, j) = (err[i].cwiseProduct(err[j])).sum();
            B(j, i) = B(i, j);  //symmetry
        }
    }

    for (int i = 0; i < B.rows() - 1; i++){
        B(B.rows() - 1, i) = -1.0;
        B(i, B.cols() - 1) = -1.0;
    }

    B(B.rows() - 1, B.cols() - 1) = 0.0;
}

double DIIS::calc_rms(Matrix &mat) const{
    double rms = 0;
    
    for (int i = 0; i < mat.rows(); i++){
        for (int j = 0; j < mat.cols(); j++){
            rms += pow(mat(i, j), 2);
        }
    }
    
    rms /= (mat.rows() * mat.cols());
    
    return sqrt(rms);
}
