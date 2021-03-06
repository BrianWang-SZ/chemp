#include "Helper.hpp"

namespace Helper{
    void print_matrix(double **mat, int size){

        for (int i = 0; i < size; i++){
            printf("%12d", i + 1);
        }

        printf("\n\n");

        for (int i = 0; i < size; i++){
            printf("%5d", i + 1);
            for (int j = 0; j < size; j++){
                printf("%12.7f", mat[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    void print_matrix(Matrix mat){
        for (int i = 0; i < mat.cols(); i++){
            printf("%12d", i + 1);
        }

        printf("\n\n");

        for (int i = 0; i < mat.rows(); i++){
            printf("%5d", i + 1);
            for (int j = 0; j < mat.cols(); j++){
                printf("%12.7f", mat(i, j));
            }
            printf("\n");
        }
        printf("\n");
    }
    
    void print_verbose(Matrix mat){
        for (int i = 0; i < mat.cols(); i++){
            printf("%21d", i + 1);
        }
        
        printf("\n\n");
        
        for (int i = 0; i < mat.rows(); i++){
            printf("%12d", i + 1);
            for (int j = 0; j < mat.cols(); j++){
                printf("%21.12e", mat(i, j));
            }
            printf("\n");
        }
        printf("\n");
    }

    double** create2d(int size){
        double **vec = new double*[size];
        
        for (int i = 0; i < size; i++){
            vec[i] = new double[size];
            for (int j = 0; j < size; j++){
                vec[i][j] = 0.0;
            }
        }
        return vec;
    }

    double calc_rms(Matrix mat, Matrix newmat){
        if (mat.rows() != newmat.rows() || mat.cols() != newmat.cols()) exit(-1);
        
        double sum = 0.0;
        for (int i = 0; i < mat.rows(); i++){
            for (int j = 0; j < mat.cols(); j++){
                sum += pow(mat(i, j) - newmat(i, j), 2);
            }
        }
        return sqrt(sum);
    }

    void free2d(double **vec, int size){
        for (int i = 0; i < size; i++){
            delete[] vec[i];
        }
        delete[] vec;
    }

    double**** create4d(int size){
        double ****mat = new double***[size];
        for (int i = 0; i < size; i++){
            mat[i] = new double**[size];
            for (int j = 0; j < size; j++){
                mat[i][j] = new double*[size];
                for (int k = 0; k < size; k++){
                    mat[i][j][k] = new double[size];
                    for (int l = 0; l < size; l++){
                        mat[i][j][k][l] = 0.0;
                    }
                }
            } 
        }
        return mat;
    }

    void free4d(double ****mat, int size){
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                for (int k = 0; k < size; k++){
                    delete[] mat[i][j][k];
                }
                delete[] mat[i][j];
            }
            delete[] mat[i];
        }
        delete[] mat;
    }
}
