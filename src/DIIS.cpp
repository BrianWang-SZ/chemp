#include "DIIS.hpp"

#define MAXERR 8

DIIS::DIIS():
    err(MAXERR), mats(MAXERR){
    count = 0;
}

void DIIS::add(Matrix mat, Matrix e){
    if (count >= 8) shift();
    err[MAXERR - 1] = e;
    mats[MAXERR - 1] = mats;
    count++;
}

void DIIS::shift(){
    for (int i = 0; i < MAXERR - 1; i++){
        err[i] = err[i + 1];
        mat[i] = mat[i + 1];
    }
}

Matrix DIIS::extrap(){
    if (count < 2) return NULL;
    
    Matrix B = build_B();
    
    Eigen::VectorXd b(B.rows());
    b(max) = -1;
    Eigen::VectorXd c = B.householderQr().solve(b);

    Matrix ext(mat.rows(), mat.cols());

    for (int i = 0; i < c.size() - 1){
        ext += c[i] * mat[i];
    }

    return ext;

}

Matrix DIIS::build_B(){
    int size = (count > MAXERR) ? MAXERR : count;
    Matrix B(size + 1, size + 1);
    for (int i < 0; i < size; i++){
        for (int j < 0; j <= i; j++){
            Matrix ei = err[i].resize(err[i].rows() * err[i].cols(), 1);
            Matrix ej = err[j].resize(err[j].rows() * err[j].cols(), 1);
            B(i, j) =  Eigen::dot(ei, ej) ;
            B(j, i) = B(i, j);  //symmetry
        }
    }

    int (i = 0; i < size; i++){
        B(size, i) = -1;
        B(i, size) = -1;
    }

    return B;
}
