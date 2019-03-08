#ifndef _DYNAMIC_2D_ARRAY_H_
#define _DYNAMIC_2D_ARRAY_H_

#include <stdio.h>
#include <stdlib.h>
             
//typedef double integer;
typedef hsize_t integer;
//typedef short sint;
typedef unsigned short sint;
typedef char cint;
typedef unsigned long tint;
typedef float real;
//typedef double real;
#define MPI_REALNUM MPI_FLOAT

integer **allocate_dynamic_2d_array_integer(integer nrows, integer ncols);
void free_dynamic_2d_array_integer(integer** array_dynamic);
void print_matrix_integer(integer** array_dynamic, integer nrows, integer ncols, char* fmt_string);

sint **allocate_dynamic_2d_array_sint(integer nrows, integer ncols);
void free_dynamic_2d_array_sint(sint** array_dynamic);
void print_matrix_sint(sint** array_dynamic, integer nrows, integer ncols, char* fmt_string);

cint **allocate_dynamic_2d_array_cint(integer nrows, integer ncols);
void free_dynamic_2d_array_cint(cint** array_dynamic);
void print_matrix_cint(cint** array_dynamic, integer nrows, integer ncols, char* fmt_string);

real **allocate_dynamic_2d_array_real(integer nrows, integer ncols);
void free_dynamic_2d_array_real(real** array_dynamic);
void print_matrix_real(real** array_dynamic, integer nrows, integer ncols, char* fmt_string);

integer** allocate_dynamic_2d_array_integer(integer nrows, integer ncols) {
    /* here is the method to correct the non-contiguous memory problem */
    integer i;
    integer** array_dynamic=(integer**)malloc(nrows*sizeof(integer*));
    integer* data=(integer*)malloc(nrows*ncols*sizeof(integer));
    for (i=0; i<nrows; i++){
        array_dynamic[i]=&(data[ncols*i]);
    }
    return array_dynamic;
}

void free_dynamic_2d_array_integer(integer** array_dynamic){
    free((void*)array_dynamic[0]);
    free((void*)array_dynamic);
}

void print_matrix_integer(integer** array_dynamic, integer nrows, integer ncols, char* fmt_string) {
    integer i,j;
    printf("\n");
    for (i = 0; i < nrows; i++){
        for (j = 0; j < ncols; j++){
            printf(fmt_string, array_dynamic[i][j]);
        }
        printf("\n");
    }
}

sint** allocate_dynamic_2d_array_sint(integer nrows, integer ncols) {
    /* here is the method to correct the non-contiguous memory problem */
    integer i;
    sint** array_dynamic=(sint**)malloc(nrows*sizeof(sint*));
    sint* data=(sint*)malloc(nrows*ncols*sizeof(sint));
    for (i=0; i<nrows; i++){
        array_dynamic[i]=&(data[ncols*i]);
    }
    return array_dynamic;
}

void free_dynamic_2d_array_sint(sint** array_dynamic){
    free((void*)array_dynamic[0]);
    free((void*)array_dynamic);
}

void print_matrix_sint(sint** array_dynamic, integer nrows, integer ncols, char* fmt_string) {
    integer i,j;
    printf("\n");
    for (i = 0; i < nrows; i++){
        for (j = 0; j < ncols; j++){
            printf(fmt_string, array_dynamic[i][j]);
        }
        printf("\n");
    }
}

cint** allocate_dynamic_2d_array_cint(integer nrows, integer ncols) {
    /* here is the method to correct the non-contiguous memory problem */
    integer i;
    cint** array_dynamic=(cint**)malloc(nrows*sizeof(cint*));
    cint* data=(cint*)malloc(nrows*ncols*sizeof(cint));
    for (i=0; i<nrows; i++){
        array_dynamic[i]=&(data[ncols*i]);
    }
    return array_dynamic;
}

void free_dynamic_2d_array_cint(cint** array_dynamic){
    free((void*)array_dynamic[0]);
    free((void*)array_dynamic);
}

void print_matrix_cint(cint** array_dynamic, integer nrows, integer ncols, char* fmt_string) {
    integer i,j;
    printf("\n");
    for (i = 0; i < nrows; i++){
        for (j = 0; j < ncols; j++){
            printf(fmt_string, array_dynamic[i][j]);
        }
        printf("\n");
    }
}

real** allocate_dynamic_2d_array_real(integer nrows, integer ncols) {
    /* here is the method to correct the non-contiguous memory problem */
    int i;
    real** array_dynamic=(real**)malloc(nrows*sizeof(real*));
    real* data=(real*)malloc(nrows*ncols*sizeof(real));
    for (i=0; i<nrows; i++){
        array_dynamic[i]=&(data[ncols*i]);
    }
    return array_dynamic;
}

void free_dynamic_2d_array_real(real** array_dynamic){
    free((void*)array_dynamic[0]);
    free((void*)array_dynamic);
}

void print_matrix_real(real** array_dynamic, integer nrows, integer ncols, char* fmt_string) {
    integer i,j;
    printf("\n");
    for (i = 0; i < nrows; i++){
        for (j = 0; j < ncols; j++){
            printf(fmt_string, array_dynamic[i][j]);
        }
        printf("\n");
    }
}

void print_matrix_1d_real(real* array_dynamic, integer dim, char* fmt_string) {
    integer i;
    printf("\n");
    for (i = 0; i < dim; i++){
        printf(fmt_string, array_dynamic[i]);
    }
    printf("\n");
}

#endif // #ifndef _DYNAMIC_2D_ARRAY_H_
