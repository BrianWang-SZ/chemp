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
}