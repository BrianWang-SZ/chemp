#include "DIIS.hpp"

#define MAXERR 8

DIIS::DIIS():
    err(MAXERR), mats(MAXERR){
    count = 0;
}

void DIIS::add(Matrix mat, Matrix e){
    if (count >= 8) shift();
    
    err[MAXERR - 1] = e;
    mats[MAXERR - 1] = mat;
    
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
    
    Eigen::VectorXd b(B.rows());
    b[B.rows() - 1] = -1;
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
            VectorXd vi = Map<VectorXd> (err[i].data(), pow(err[i].size(), 2));
            VectorXd vj = Map<VectorXd> (err[j].data(), pow(err[j].size(), 2));
            B(i, j) = vi.dot(vj);
            B(j, i) = B(i, j);  //symmetry
        }
    }

    for (int i = 0; i < size; i++){
        B(size, i) = -1;
        B(i, size) = -1;
    }

    return B;
}
