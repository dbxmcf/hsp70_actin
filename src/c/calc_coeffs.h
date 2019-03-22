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
                     real *sum_min, real *sum_max,
                     real *sum_min_jac, real *sum_max_jac, real* num_sim)
{
    tint i;
    sint c_min,c_max,cj_min,cj_max,num_c_sim;  
    real sum_c_min=0.0,sum_c_max=0.0,sum_cj_min=0.0,sum_cj_max=0.0,sum_num_sim=0.0;
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
    *sum_min = sum_c_min;
    *sum_max = sum_c_max;
    *sum_min_jac = sum_cj_min;
    *sum_max_jac = sum_cj_max;
    *num_sim = sum_num_sim;
    //printf("%7.3f\n",sum);
    //return sum;
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
    tint i, j, idx_a, idx_b;
    sint *a, *b;
    //cint *a_jac, *b_jac;
    real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
    real numerator_sarika, denomenator_sarika, dist_sarika;
    real num_sim, numerator_jac, denomenator_jac, numerator_gen_jac, denomenator_gen_jac;
    real a_sum, b_sum, one_an, one_bn, result;

    // perform full permutation of all pairs of part a and part b
    if (part_a_dim1 != part_b_dim1) {
        printf("dim1s do not match!");
        exit(1);
    }

    tint dim1 = part_a_dim1;
    // calculate some preparation values
    real *restrict data_sum_a = (real*)malloc(part_a_dim0*sizeof(real));
    real *restrict data_sum_b = (real*)malloc(part_b_dim0*sizeof(real));
    real *restrict one_data_norm_a = (real*)malloc(part_a_dim0*sizeof(real));
    real *restrict one_data_norm_b = (real*)malloc(part_b_dim0*sizeof(real));

    //tint **restrict data_jac_a = allocate_dynamic_2d_array_integer(part_a_dim0, dim1);
    //tint **restrict data_jac_b = allocate_dynamic_2d_array_integer(part_b_dim0, dim1);
    //cint **restrict data_jac_a = allocate_dynamic_2d_array_cint(part_a_dim0, dim1);
    //cint **restrict data_jac_b = allocate_dynamic_2d_array_cint(part_b_dim0, dim1);

    //real *restrict normal = (real*)malloc(part_a_dim0*part_b_dim0*sizeof(real));
    //real *restrict generalised = (real*)malloc(part_a_dim0*part_b_dim0*sizeof(real));
    //real *restrict sarika = (real*)malloc(part_a_dim0*part_b_dim0*sizeof(real));
    //real *restrict wu = (real*)malloc(part_a_dim0*part_b_dim0*sizeof(real));
    //real *restrict cosine = (real*)malloc(part_a_dim0*part_b_dim0*sizeof(real));

    //tint **cmb_ab = allocate_dynamic_2d_array_integer(part_a_dim0, part_b_dim0);
    //normal[idx_a,idx_b] = dist_jac
    //normal[idx_b,idx_a] = dist_jac
    //generalised[idx_a,idx_b] = dist_gen_jac
    //generalised[idx_b,idx_a] = dist_gen_jac
    //sarika[idx_a,idx_b] = dist_sarika
    //sarika[idx_b,idx_a] = dist_sarika
    //wu[idx_a,idx_b] = dist_wu
    //wu[idx_b,idx_a] = dist_wu
    //cosine[idx_a,idx_b] = result*100
    //cosine[idx_b,idx_a] = result*100

    // part_a values
    //#pragma acc kernels
    for (i=0;i<part_a_dim0;i++) {
        data_sum_a[i] = 0;
        for (j=0;j<part_a_dim1;j++) {
            //data_jac_a[i][j]=0;
            //if (data_part_a[i][j]>0) 
            //    data_jac_a[i][j]=1;
            data_sum_a[i] += data_part_a[i][j];
        }
        //one_data_norm_a[i] = 1.0/vec_norm(data_part_a[i], dim1);
        one_data_norm_a[i] = 0.0;
    }

    // part_b values
    //#pragma acc kernels
    for (i=0;i<part_b_dim0;i++) {
        data_sum_b[i] = 0;
        for (j=0;j<part_b_dim1;j++) {
            //data_jac_b[i][j]=0;
            //if (data_part_b[i][j]>0) 
            //    data_jac_b[i][j]=1;
            data_sum_b[i] += data_part_b[i][j];
        }
        //one_data_norm_b[i] = 1.0/vec_norm(data_part_b[i], dim1);
        one_data_norm_b[i] = 0.0;
    }

    /* divide the large cpu memory data block to small device blocks

     */
    tint device_block_id, device_block_num;
    tint device_num_data_chunks_parts = 3;
    //tint device_num_data_chunks_part_b = 3;
    integer **device_cmbs=NULL;
    integer device_num_cmbs = 0;
    tint device_blk_part_a_avg_lines = part_a_dim0/device_num_data_chunks_parts;
    tint device_blk_part_b_avg_lines = part_b_dim0/device_num_data_chunks_parts;
    device_block_num = device_num_data_chunks_parts*device_num_data_chunks_parts;

    
    //device_cmbs=combination_util(device_num_data_chunks,&device_num_cmbs); 
    //printf("num_cmb=%d,num_cmbs1=%d\n",num_cmbs,num_cmbs1);
    //print_matrix(cmbs, num_cmbs, 2, "%3d");
    integer device_chunk0, device_chunk1;

    if (mpi_rank < num_cmbs) { // an off-diagnal full block
        //printf("mpi_rank=%d\n",mpi_rank);
        mpi_rk_chunk0 = cmbs[mpi_rank][0];
        mpi_rk_chunk1 = cmbs[mpi_rank][1];
    }
    else { // two diagnal triangles
        mpi_rk_chunk0 = (mpi_rank-num_cmbs)*2;
        mpi_rk_chunk1 = mpi_rk_chunk0 + 1;
    }
    if (verbose)
        printf("mrk[%d]:mpi_rk_chunk0=%d, mpi_rk_chunk1=%d\n",mpi_rank,mpi_rk_chunk0, mpi_rk_chunk1);
    free_dynamic_2d_array_integer(cmbs);


    tint idx_out = 0;

    tint device_blk_start, device_blk_end;
    for (device_block_id=0;device_block_id<device_block_num;device_block_id++)
    {
        device_blk_start_part_a = device_blk_start_array_a[device_block_id];
        device_blk_start_part_b = device_blk_start_array_b[device_block_id];
        device_blk_size_part_a = device_blk_size_array_a[device_block_id];
        device_blk_size_part_b = device_blk_size_array_b[device_block_id];
    // the large loop that calculates the matrix
//#pragma acc data \
//    copy(data_part_a[0:part_a_dim0][0:part_a_dim1],\
//            data_part_b[0:part_b_dim0][0:part_b_dim1],\
//            data_sum_a[0:part_a_dim0],\
//            data_sum_b[0:part_b_dim0],\
//            one_data_norm_a[0:part_a_dim0],\
//            one_data_norm_b[0:part_b_dim0])
#pragma acc data \
    copy(device_data_part_a[0:device_blk_size_a][0:part_a_dim1],\
            device_data_part_b[0:device_blk_size_b][0:part_b_dim1],\
            device_data_sum_a[0:device_blk_size_a],\
            device_data_sum_b[0:device_blk_size_b],\
            device_one_data_norm_a[0:device_blk_size_a],\
            device_one_data_norm_b[0:device_blk_size_b])
    {

        for (idx_a=0;idx_a<device_blk_size_a;idx_a++) {
            for (idx_b=0;idx_b<device_blk_size_b;idx_b++){
                // 
                a = device_data_part_a[idx_a];
                a_sum = device_data_sum_a[idx_a];
                b = device_data_part_b[idx_b];
                b_sum = device_data_sum_b[idx_b];

                sum_min_max_vec(a, b, dim1, 
                                &numerator_gen_jac,&denomenator_gen_jac,
                                &numerator_jac,&denomenator_jac,&num_sim);

                //one_an = one_data_norm_a[idx_a];
                //one_bn = one_data_norm_b[idx_b];
                //result = 1.0 - vec_dot(a,b,dim1)*one_an*one_bn;
                //result = vec_dot(a,b,dim1)*one_an*one_bn;

                dist_gen_jac = 1.0-numerator_gen_jac/denomenator_gen_jac;

                dist_jac = 1.0-numerator_jac/denomenator_jac;

                denomenator_wu = MIN(denomenator_gen_jac,MAX(a_sum,b_sum));
                dist_wu = 1.0-numerator_gen_jac/denomenator_wu;

                numerator_sarika = num_sim;
                denomenator_sarika = a_sum+b_sum;
                dist_sarika = 1.0-numerator_sarika/denomenator_sarika;

                rp->normal[idx_out] = dist_jac;
                rp->generalised[idx_out] = dist_gen_jac;
                rp->sarika[idx_out] = dist_sarika;
                rp->wu[idx_out] = dist_wu;
                rp->cosine[idx_out] = 0.0; //result*100;
                idx_out++;
            }
            }
        }
        } // for device_block_id

        free(data_sum_a);
        free(data_sum_b);
        free(one_data_norm_a);
        free(one_data_norm_b);

        return 0;
    }



    int calc_coeffs_diagnol_triangle(sint **restrict data, tint dim0, tint dim1,
            result_pointers_diagnol *rpd)
    {
        tint i, j, idx_a, idx_b;
        sint *a, *b;
        real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
        real numerator_sarika, denomenator_sarika, dist_sarika;
        real num_sim, numerator_jac, denomenator_jac, numerator_gen_jac, denomenator_gen_jac;
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
        copy(data[0:dim0][0:dim1],\
                data_sum[0:dim0],\
                one_data_norm[0:dim0])
        {

            for (i=0;i< num_cmbs;i++) {
                idx_a = cmbs[i][0], idx_b = cmbs[i][1];

                a = data[idx_a];
                a_sum = data_sum[idx_a];
                b = data[idx_b];
                b_sum = data_sum[idx_b];

                sum_min_max_vec(a, b, dim1, 
                                &numerator_gen_jac,&denomenator_gen_jac,
                                &numerator_jac,&denomenator_jac,&num_sim);

                //one_an = one_data_norm[idx_a];
                //one_bn = one_data_norm[idx_b];
                //real adotb = vec_dot(a,b,dim1);
                //printf("%d %d %.3e\n", idx_a, idx_b, adotb);
                //result = 1.0 - adotb*one_an*one_bn;
                //result = adotb*one_an*one_bn;

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
                rpd->cosine[idx_rpd] = 0.0; //result*100;
                //if (result >100 || result <0) {
                //    printf("adotb=%7.3e,one_an=%7.3e,one_bn=%7.3e\n",adotb,one_an,one_bn);
                //}
                idx_rpd++;
            }
        }

        /* pass out the data */
        //for (i=0;i<idx_rpd;i++) {
        //    //printf("%7.3f ",normal[i]);
        //    rpd->normal[i] = normal[i];
        //    rpd->generalised[i] = generalised[i];
        //    rpd->sarika[i] = sarika[i];
        //    rpd->wu[i] = wu[i];
        //    rpd->cosine[i] = cosine[i];
        //}

        //free(normal);
        //free(generalised);
        //free(sarika);
        //free(wu);
        //free(cosine);

        free_dynamic_2d_array_integer(cmbs);
        //free_dynamic_2d_array_cint(data_jac);
        free(data_sum);
        free(one_data_norm);

        return 0;

    }

#endif
