#ifndef _CALC_COEFFS_H_
#define _CALC_COEFFS_H_
#include <stdio.h>
#include <math.h>
#include <omp.h>

//#define N 3
//#define D 5
//#define C 3

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


//typedef float real;
//typedef unsigned short tint;

typedef struct result_arrays {
    real **normal;
    real **generalised;
    real **wu;
    real **sarika;
    real **cosine;
} result_pointers;

typedef struct result_arrays_diagnol {
    real *normal;
    real *generalised;
    real *wu;
    real *sarika;
    real *cosine;
    char *normal_name;
    char *generalised_name;
    char *wu_name;
    char *sarika_name;
    char *cosine_name;
    integer vec_dim;
    integer total_lines;
    integer start_loc;
    integer chunk_start_a;
    integer chunk_start_b;
    integer chunk_count_a;
    integer chunk_count_b;

} result_pointers_diagnol;

void sum_min_max_vec(sint *restrict a, sint *restrict b, tint vec_dim, 
        real a_sum, real b_sum, result_pointers_diagnol *rpd, tint idx_rpd)
//        real *sum_min, real *sum_max,
//        real *sum_min_jac, real *sum_max_jac, real* num_sim)
{
    tint i;
    sint c_min,c_max,cj_min,cj_max,num_c_sim;  
    real sum_c_min=0.0,sum_c_max=0.0,sum_cj_min=0.0,sum_cj_max=0.0,sum_num_sim=0.0;
    real numerator_gen_jac, denomenator_gen_jac, numerator_jac, denomenator_jac, num_sim;
    real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
    real numerator_sarika, denomenator_sarika, dist_sarika;
    real one_an, one_bn, result;

#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(c_min,c_max,cj_min,cj_max,i) \
    reduction(+:sum_c_min,sum_c_max,sum_cj_min,sum_cj_max,sum_num_sim)
    for (i=0;i<vec_dim;i++) {
        // c = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // min(x, y)
        c_min = ((a[i])<(b[i]))?(a[i]):(b[i]);
        c_max = ((a[i])>(b[i]))?(a[i]):(b[i]);
        cj_min = (c_min>0)?1:0;
        cj_max = (c_max>0)?1:0;
        //printf("%d, %d, %d\n", a[i],b[i],c);
        //if (cj_min>0) num_c_sim = a[i]+b[i];
        num_c_sim = (c_min>0)?(a[i]+b[i]):0;
        sum_c_min += c_min;
        sum_c_max += c_max;
        sum_cj_min += cj_min;
        sum_cj_max += cj_max;
        sum_num_sim += num_c_sim;
    }
    //   *sum_min = sum_c_min;
    //   *sum_max = sum_c_max;
    //   *sum_min_jac = sum_cj_min;
    //   *sum_max_jac = sum_cj_max;
    //   *num_sim = sum_num_sim;

    //void sum_min_max_vec(sint *restrict a, sint *restrict b, tint vec_dim, 
    //        real *sum_min, real *sum_max,
    //        real *sum_min_jac, real *sum_max_jac, real* num_sim)

    //sum_min_max_vec(a, b, dim1, 
    //                    &numerator_gen_jac,&denomenator_gen_jac,
    //                    &numerator_jac,&denomenator_jac,&num_sim);

    numerator_gen_jac = sum_c_min;
    denomenator_gen_jac = sum_c_max;
    numerator_jac = sum_cj_min;
    denomenator_jac = sum_cj_max;
    num_sim = sum_num_sim;

    result = 1.0; //adotb*one_an*one_bn;
    dist_gen_jac = 1.0-numerator_gen_jac/denomenator_gen_jac;
    dist_jac = 1.0-numerator_jac/denomenator_jac;

    denomenator_wu = MIN(denomenator_gen_jac,MAX(a_sum,b_sum));
    dist_wu = 1.0-numerator_gen_jac/denomenator_wu;

    numerator_sarika = num_sim;
    denomenator_sarika = a_sum+b_sum;
    dist_sarika = 1.0-numerator_sarika/denomenator_sarika;

    rpd->normal[idx_rpd] = dist_jac;
    rpd->generalised[idx_rpd] = dist_gen_jac;
    rpd->sarika[idx_rpd] = dist_sarika;
    rpd->wu[idx_rpd] = dist_wu;
    rpd->cosine[idx_rpd] = result*100;
}

// https://www.geeksforgeeks.org/compute-the-minimum-or-maximum-max-of-two-integers-without-branching/
// https://stackoverflow.com/questions/24529504/find-out-max-min-of-two-number-without-using-if-else
// http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
real sum_minimum_vec(sint *restrict a, sint *restrict b, tint vec_dim, real *sum_jac)
{
    tint i;
    sint c,cj;  
    real sum=0,sum_cj=0;
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(c,i) reduction(+:sum,sum_cj)
    for (i=0;i<vec_dim;i++) {
        // c = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // min(x, y)
        c = ((a[i])<(b[i]))?(a[i]):(b[i]);
        cj = (c>0)?1:0;
        //printf("%d, %d, %d\n", a[i],b[i],c);
        sum += c;
        sum_cj +=cj;
    }
    *sum_jac = sum_cj;
    //printf("%7.3f\n",sum);
    return sum;
}

real sum_maximum_vec(sint *restrict a, sint *restrict b, tint vec_dim, real *sum_jac)
{
    tint i;
    sint c,cj;  
    real sum=0,sum_cj=0;  
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(c,cj,i) reduction(+:sum,sum_cj)
    for (i=0;i<vec_dim;i++) {
        //c = a[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // max(x, y)
        c = ((a[i])>(b[i]))?(a[i]):(b[i]);
        cj = (c>0)?1:0;
        sum += c;
        sum_cj +=cj;
    }
    *sum_jac = sum_cj;
    return sum;
}

real sum_minimum_vec_jac(sint *restrict a, sint *restrict b, tint vec_dim)
{
    tint i;
    sint c; //,c1;  
    real sum=0;
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(c,i) reduction(+:sum)
    for (i=0;i<vec_dim;i++) {
        // c = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // min(x, y)
        //c = (a[i] && b[i]);
        //c1 = ((a[i])<(b[i]))?(a[i]):(b[i]);
        //if (c1>0) c1=1;

        c = ((a[i])<(b[i]))?(a[i]):(b[i]);
        if (c>0) c=1;

        //c=1;
        //if (!(a[i]>0)) c=0;
        //if (!(b[i]>0)) c=0;
        //c = (a[i]>0 && b[i]>0)?1:0;
        //if (c1 != c) {
        //    printf("a=%d, b=%d, c1=%d, c=%d\n", a[i],b[i],c1,c);
        //    exit(0);
        //}
        sum += c;
    }
    //printf("%7.3f\n",sum);
    return sum;
}

real sum_maximum_vec_jac(sint *restrict a, sint *restrict b, tint vec_dim)
{
    tint i;
    sint c;  
    real sum=0;  
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(c,i) reduction(+:sum)
    for (i=0;i<vec_dim;i++) {
        //c = a[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // max(x, y)
        //c = (a[i] || b[i]);
        //c = (a[i]>0 || b[i]>0)?1:0;

        c = ((a[i])>(b[i]))?(a[i]):(b[i]);
        if (c>0) c=1;

        //c=0;
        //if (a[i]>0) c=1;
        //if (b[i]>0) c=1;
        sum += c;
    }
    return sum;
}

real sum_minimum_vec_cint(cint *restrict a, cint *restrict b, tint vec_dim)
{
    tint i, c;  
    real sum=0;
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(c,i) reduction(+:sum)
    for (i=0;i<vec_dim;i++) {
        // c = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // min(x, y)
        c = ((a[i])<(b[i]))?(a[i]):(b[i]);
        //printf("%d, %d, %d\n", a[i],b[i],c);
        sum += c;
    }
    //printf("%7.3f\n",sum);
    return sum;
}

real sum_maximum_vec_cint(cint *restrict a, cint *restrict b, tint vec_dim)
{
    tint i, c;  
    real sum=0;  
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(c,i) reduction(+:sum)
    for (i=0;i<vec_dim;i++) {
        //c = a[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // max(x, y)
        c = ((a[i])>(b[i]))?(a[i]):(b[i]);
        sum += c;
    }
    return sum;
}

real get_non_zeros_pair(sint *restrict a, sint *restrict b, tint vec_dim)
{
    real sum = 0;
    tint i;
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(i) reduction(+:sum)
    for (i=0;i<vec_dim;i++) {
        if (a[i]>0 && b[i]>0)
            //if (a[i] & b[i])
            sum += a[i]+b[i];
    }
    return sum;
}

tint get_non_zeros_single(tint *a, tint vec_dim)
{
    tint count = 0;
    tint i;
    for (i=0;i<vec_dim;i++)
        if (a[i] != 0)
            count++;
    return count;
}

void get_sum_vec(tint *a, tint *b, tint vec_dim)
{
    real sum=0;
    tint i;
    for (i=0;i<vec_dim;i++)
        sum += a[i]+b[i];
}

void vec_add(tint *a, tint *b, tint *c,tint vec_dim)
{
    tint i;  
    for(i=0; i<vec_dim; i++) 
        c[i] = a[i]+b[i];
}

void vec_and(tint *a, tint *b, tint *c,tint vec_dim)
{
    tint i;  
    for(i=0; i<vec_dim; i++) 
        c[i] = a[i] & b[i];
}

void vec_or(tint *a, tint *b, tint *c,tint vec_dim)
{
    tint i;  
    for(i=0; i<vec_dim; i++) 
        c[i] = a[i] | b[i];
}

real vec_norm(sint *restrict a, sint vec_dim)
{
    tint i;
    real norm=0;
#pragma acc parallel loop present(a[0:vec_dim])
#pragma omp parallel for private(i) reduction(+:norm)
    for(i=0;i<vec_dim; i++) {
        norm += a[i]*a[i];
    }
    return sqrt(norm);
}

real vec_dot(sint *restrict a, sint *restrict b, tint vec_dim)
{
    tint i;
    real sum=0; 
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(i) reduction(+:sum)
    for(i=0;i<vec_dim; i++) {
        sum += a[i]*b[i];
    }
    return sum;
}

// https://stackoverflow.com/questions/20013693/read-csv-file-to-a-2d-array-on-c
void read_h5()
{

}

int calc_coeffs_off_diagnol_block(sint **restrict data_part_a, tint part_a_dim0, tint part_a_dim1,
        sint **restrict data_part_b, tint part_b_dim0, tint part_b_dim1,
        result_pointers_diagnol *rp)
{
    //int i, j, idx_a, idx_b;
    int is_diagnol = 0;
    int is_dvc_blk_diagnol = 0;
    tint idx_b_loop_start = 0;
    if (data_part_a == data_part_b) {
        is_diagnol = 1;
    }
    //printf("is_diagnol=%ld\n",is_diagnol);
    tint i, j, idx_a, idx_b;
    //sint *a, *b;
    //cint *a_jac, *b_jac;
    //real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
    //real numerator_sarika, denomenator_sarika, dist_sarika;
    //real num_sim, numerator_jac, denomenator_jac, numerator_gen_jac, denomenator_gen_jac;
    //real a_sum, b_sum, one_an, one_bn, result;

    // perform full permutation of all pairs of part a and part b
    if (part_a_dim1 != part_b_dim1) {
        printf("dim1s do not match!");
        exit(1);
    }

    tint dim1 = part_a_dim1;
    // calculate some preparation values
    real *restrict data_sum_a = (real*)malloc(part_a_dim0*sizeof(real));
    real *restrict one_data_norm_a = (real*)malloc(part_a_dim0*sizeof(real));
    real *restrict data_sum_b = NULL;
    real *restrict one_data_norm_b = NULL;

    // part_a values, for diagnol process only part_a is needed
    //#pragma acc kernels
    for (i=0;i<part_a_dim0;i++) {
        data_sum_a[i] = 0;
        for (j=0;j<part_a_dim1;j++) {
            data_sum_a[i] += data_part_a[i][j];
        }
        //one_data_norm_a[i] = 1.0/vec_norm(data_part_a[i], dim1);
        one_data_norm_a[i] = 0.0;
    }

    // part_b values
    if (is_diagnol) {
        data_sum_b = data_sum_a;
        one_data_norm_b = one_data_norm_a;
    }
    else {
        data_sum_b = (real*)malloc(part_b_dim0*sizeof(real));
        one_data_norm_b = (real*)malloc(part_b_dim0*sizeof(real));
        //#pragma acc kernels
        for (i=0;i<part_b_dim0;i++) {
            data_sum_b[i] = 0;
            for (j=0;j<part_b_dim1;j++) {
                data_sum_b[i] += data_part_b[i][j];
            }
            //one_data_norm_b[i] = 1.0/vec_norm(data_part_b[i], dim1);
            one_data_norm_b[i] = 0.0;
        }
    }

    tint idx_out = 0;
    /* can they be different*/
    tint dvc_blk_part_a_num = device_block_parts_num, dvc_blk_part_b_num = dvc_blk_part_a_num;
    if (debug_info)
        printf("number of device blocks=%ld\n",dvc_blk_part_a_num);

    tint *dvc_blk_part_a_start = NULL, *dvc_blk_part_a_size  = NULL;
    tint *dvc_blk_part_b_start = NULL, *dvc_blk_part_b_size  = NULL;

    distribute_parts_start_size(part_a_dim0, dvc_blk_part_a_num, &dvc_blk_part_a_start, &dvc_blk_part_a_size);
    distribute_parts_start_size(part_b_dim0, dvc_blk_part_b_num, &dvc_blk_part_b_start, &dvc_blk_part_b_size);

    sint **restrict dvc_blk_part_a;
    sint **restrict dvc_blk_part_b;
    tint dvc_blk_part_a_start_idx, dvc_blk_part_b_start_idx;
    tint dvc_blk_part_b_loop_begin;
    tint dvc_blk_part_a_dim0, dvc_blk_part_b_dim0;
    tint global_idx_a,global_idx_b;
    real *restrict dvc_blk_sum_a; // = (real*)malloc(part_a_dim0*sizeof(real));
    real *restrict dvc_blk_sum_b;

    tint idx_dvc_blk_part_a, idx_dvc_blk_part_b;

    //return 0;

    //dvc_blk_part_b_loop_begin = 0;
    /* outer loop for the device blocks*/
    for (idx_dvc_blk_part_a=0;idx_dvc_blk_part_a<dvc_blk_part_a_num;idx_dvc_blk_part_a++) {
        dvc_blk_part_b_loop_begin = 0;
        if (is_diagnol) {
            dvc_blk_part_b_loop_begin = idx_dvc_blk_part_a;
        }
        for (idx_dvc_blk_part_b=dvc_blk_part_b_loop_begin;
             idx_dvc_blk_part_b<dvc_blk_part_b_num;
             idx_dvc_blk_part_b++) {
                            
            /* get the device block start and size, part a */
            dvc_blk_part_a_start_idx = dvc_blk_part_a_start[idx_dvc_blk_part_a]; /* get the index */
            dvc_blk_part_a = &data_part_a[dvc_blk_part_a_start_idx]; /* get the address */
            dvc_blk_part_a_dim0 = dvc_blk_part_a_size[idx_dvc_blk_part_a];
            //printf("dvc_blk_part_a_start_idx-%ld\n",dvc_blk_part_a_start_idx);
            //printf("dvc_blk_part_a_dim0-%ld\n",dvc_blk_part_a_dim0);
            /* part b */
            dvc_blk_part_b_start_idx = dvc_blk_part_b_start[idx_dvc_blk_part_b];
            dvc_blk_part_b = &data_part_b[dvc_blk_part_b_start_idx];
            dvc_blk_part_b_dim0 = dvc_blk_part_b_size[idx_dvc_blk_part_b];
            //printf("dvc_blk_part_b_start_idx-%ld\n",dvc_blk_part_b_start_idx);
            //printf("dvc_blk_part_b_dim0-%ld\n",dvc_blk_part_b_dim0);

            dvc_blk_sum_a = &data_sum_a[dvc_blk_part_a_start[idx_dvc_blk_part_a]];
            dvc_blk_sum_b = &data_sum_b[dvc_blk_part_b_start[idx_dvc_blk_part_b]];

            //break;

            is_dvc_blk_diagnol = 0;
            if ( is_diagnol && idx_dvc_blk_part_a == idx_dvc_blk_part_b )
                is_dvc_blk_diagnol = 1;

                /* within block loop */
                #pragma acc data \
                copyin(dvc_blk_part_a[0:dvc_blk_part_a_dim0][0:part_a_dim1],\
                       dvc_blk_part_b[0:dvc_blk_part_b_dim0][0:part_b_dim1])
                {

                    for (idx_a=0;idx_a<dvc_blk_part_a_dim0;idx_a++) {
                        idx_b_loop_start =0;
                        if (is_dvc_blk_diagnol) {
                            idx_b_loop_start = idx_a + 1;
                        }
                        for (idx_b=idx_b_loop_start;idx_b<dvc_blk_part_b_dim0;idx_b++){
//
                            global_idx_a = dvc_blk_part_a_start_idx + idx_a;
                            global_idx_b = dvc_blk_part_b_start_idx + idx_b;

                            idx_out = global_idx_a*part_b_dim0 + global_idx_b;
                            if (is_diagnol) {
                                idx_out -= (global_idx_a+1)*(global_idx_a+2)/2;
                            }
//
                            sum_min_max_vec(dvc_blk_part_a[idx_a], dvc_blk_part_b[idx_b], dim1, 
                                            dvc_blk_sum_a[idx_a], dvc_blk_sum_b[idx_b], 
                                            rp, idx_out);
                            //printf("idx_out_tri=%ld\n",idx_out);
                            //idx_out++;
                        }
                    }
                }

        }
    }


    // the large loop that calculates the matrix
//#pragma acc data \
//    copy(data_part_a[0:part_a_dim0][0:part_a_dim1],\
//         data_part_b[0:part_b_dim0][0:part_b_dim1])
//    {
//
//        for (idx_a=0;idx_a<part_a_dim0;idx_a++) {
//            for (idx_b=0;idx_b<part_b_dim0;idx_b++){
//                // 
//                a = data_part_a[idx_a];
//                a_sum = data_sum_a[idx_a];
//                b = data_part_b[idx_b];
//                b_sum = data_sum_b[idx_b];
//
//                sum_min_max_vec(a, b, dim1, a_sum, b_sum, rp, idx_out);
//
//                idx_out++;
//            }
//        }
//    }

    free(dvc_blk_part_a_start);
    free(dvc_blk_part_a_size);
    free(dvc_blk_part_b_start);
    free(dvc_blk_part_b_size);

    free(data_sum_a);
    free(one_data_norm_a);
    if ( !is_diagnol ) { /* only off-diagnol blocks need to free part_b*/
        free(data_sum_b);
        free(one_data_norm_b);
    }

    return 0;
}



int calc_coeffs_diagnol_triangle(sint **restrict data, tint dim0, tint dim1,
        result_pointers_diagnol *rpd)
{
    tint i, j, idx_a, idx_b;
    sint *a, *b;
    //cint *a_jac, *b_jac;

    real a_sum, b_sum, one_an, one_bn, result;

    // calculate some preparation values
    real *restrict data_sum = (real*)malloc(dim0*sizeof(real));
    real *restrict one_data_norm = (real*)malloc(dim0*sizeof(real));
    integer **restrict cmbs=NULL;

    integer num_cmbs = 0, idx_rpd = 0;
    cmbs = combination_util(dim0,&num_cmbs); 
    // preparation values
    //#pragma acc kernels
    for (i=0;i<dim0;i++) {
        data_sum[i] = 0;
        for (j=0;j<dim1;j++) {
            data_sum[i] += data[i][j];
        }
        //one_data_norm[i] = 1.0/vec_norm(data[i], dim1);
    }

#pragma acc data \
    copy(data[0:dim0][0:dim1])
    {

        for (i=0;i< num_cmbs;i++) {
            idx_a = cmbs[i][0], idx_b = cmbs[i][1];

            a = data[idx_a];
            a_sum = data_sum[idx_a];
            b = data[idx_b];
            b_sum = data_sum[idx_b];

            sum_min_max_vec(a, b, dim1, a_sum, b_sum, rpd, idx_rpd);

            idx_rpd++;
        }
    }

    printf("idx_rpd=%ld\n",idx_rpd);

    /* pass out the data */
    free_dynamic_2d_array_integer(cmbs);
    free(data_sum);
    free(one_data_norm);

    return 0;

}

#endif
