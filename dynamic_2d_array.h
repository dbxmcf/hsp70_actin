#ifndef _DYNAMIC_2D_ARRAY_H_
#define _DYNAMIC_2D_ARRAY_H_

#include <stdio.h>
#include <stdlib.h>
             
//typedef double real;
typedef float real;
#define MPI_REALNUM MPI_FLOAT

real **allocate_dynamic_2d_array(int nrows, int ncols);
void free_dynamic_2d_array(real** array_dynamic);
void print_matrix(real** array_dynamic, int nrows, int ncols, char* fmt_string);

real** allocate_dynamic_2d_array(int nrows, int ncols) {
    /* here is the method to correct the non-contiguous memory problem */
    int i;
    real** array_dynamic=(real**)malloc(nrows*sizeof(real*));
    real* data=(real*)malloc(nrows*ncols*sizeof(real));
    for (i=0; i<nrows; i++){
        array_dynamic[i]=&(data[ncols*i]);
    }
    return array_dynamic;
}

void free_dynamic_2d_array(real** array_dynamic){
    free((void*)array_dynamic[0]);
    free((void*)array_dynamic);
}

void print_matrix(real** array_dynamic, int nrows, int ncols, char* fmt_string) {
    int i,j;
    printf("\n");
    for (i = 0; i < nrows; i++){
        for (j = 0; j < ncols; j++){
            printf(fmt_string, array_dynamic[i][j]);
        }
        printf("\n");
    }
}

#endif // #ifndef _DYNAMIC_2D_ARRAY_H_
