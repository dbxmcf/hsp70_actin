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
typedef int tint;

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
    tint vec_dim;
    tint total_lines;
    tint start_loc;
    tint chunk_start_a;
    tint chunk_start_b;
    tint chunk_count_a;
    tint chunk_count_b;

} result_pointers_diagnol;

// https://www.geeksforgeeks.org/compute-the-minimum-or-maximum-max-of-two-integers-without-branching/
// https://stackoverflow.com/questions/24529504/find-out-max-min-of-two-number-without-using-if-else
// http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
real sum_minimum_vec(tint *restrict a, tint *restrict b, tint vec_dim)
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

real sum_maximum_vec(tint *restrict a, tint *restrict b, tint vec_dim)
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

real get_non_zeros_pair(tint *restrict a, tint *restrict b, tint vec_dim)
{
    real sum = 0;
    tint i;
#pragma acc parallel loop present(a[0:vec_dim],b[0:vec_dim])
#pragma omp parallel for private(i) reduction(+:sum)
    for (i=0;i<vec_dim;i++) {
        if (a[i]>0 && b[i]>0)
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

real vec_norm(tint *restrict a, tint vec_dim)
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

real vec_dot(tint *restrict a, tint *restrict b, tint vec_dim)
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

int calc_coeffs_off_diagnol_block(tint **restrict data_part_a, tint part_a_dim0, tint part_a_dim1,
        tint **restrict data_part_b, tint part_b_dim0, tint part_b_dim1,
        result_pointers_diagnol *rp)
{
    //int i, j, idx_a, idx_b;
    tint i, j, idx_a, idx_b;
    tint *a, *b, *a_jac, *b_jac;
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

    tint **restrict data_jac_a = allocate_dynamic_2d_array_integer(part_a_dim0, dim1);
    tint **restrict data_jac_b = allocate_dynamic_2d_array_integer(part_b_dim0, dim1);

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
            data_jac_a[i][j]=0;
            if (data_part_a[i][j]>0) 
                data_jac_a[i][j]=1;
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
            data_jac_b[i][j]=0;
            if (data_part_b[i][j]>0) 
                data_jac_b[i][j]=1;
            data_sum_b[i] += data_part_b[i][j];
        }
        //one_data_norm_b[i] = 1.0/vec_norm(data_part_b[i], dim1);
        one_data_norm_b[i] = 0.0;
    }


    int idx_out = 0;
    // the large loop that calculates the matrix
#pragma acc data \
    copy(data_part_a[0:part_a_dim0][0:part_a_dim1],\
            data_part_b[0:part_b_dim0][0:part_b_dim1],\
            data_sum_a[0:part_a_dim0],\
            data_sum_b[0:part_b_dim0],\
            data_jac_a[0:part_a_dim0][0:part_a_dim1],\
            data_jac_b[0:part_b_dim0][0:part_b_dim1],\
            one_data_norm_a[0:part_a_dim0],\
            one_data_norm_b[0:part_b_dim0])
    {

        for (idx_a=0;idx_a<part_a_dim0;idx_a++) {
            for (idx_b=0;idx_b<part_b_dim0;idx_b++){
                // 
                a = data_part_a[idx_a];
                a_sum = data_sum_a[idx_a];
                a_jac = data_jac_a[idx_a];
                b = data_part_b[idx_b];
                b_sum = data_sum_b[idx_b];
                b_jac = data_jac_b[idx_b];

                //vec_add(a, b, summed_array, dim1);
                numerator_jac = sum_minimum_vec(a_jac, b_jac, dim1);

                denomenator_jac = sum_maximum_vec(a_jac, b_jac, dim1);
                //printf("numerator_jac=%f\n",numerator_jac);
                //printf("denomenator_jac=%f\n",denomenator_jac);
                numerator_gen_jac = sum_minimum_vec(a, b, dim1);
                denomenator_gen_jac = sum_maximum_vec(a, b, dim1);

                //printf("numerator_gen_jac=%f\n",numerator_gen_jac);
                //printf("denomenator_gen_jac=%f\n",denomenator_gen_jac);

                num_sim = get_non_zeros_pair(a, b, dim1);
                //printf("num_sim=%f\n",num_sim);

                one_an = one_data_norm_a[idx_a];
                one_bn = one_data_norm_b[idx_b];
                //result = 1.0 - vec_dot(a,b,dim1)*one_an*one_bn;
                result = vec_dot(a,b,dim1)*one_an*one_bn;

                dist_gen_jac = 1.0-numerator_gen_jac/denomenator_gen_jac;
                //if (numerator_jac < 0 || denomenator_jac > 100) {
                //if (numerator_jac < 0 ) {
                //    printf("a[0]=%d\n",a[0]);
                //    printf("b[0]=%d\n",b[0]);
                //    printf("numerator_jac:%lf \n",numerator_jac);
                //    printf("denomenator_jac:%lf \n",denomenator_jac);
                //}

                dist_jac = 1.0-numerator_jac/denomenator_jac;
                //printf("numerator_jac:%lf \n",numerator_jac);
                //printf("denomenator_jac:%lf \n",denomenator_jac);
                //printf("dist_jac:%lf \n",dist_jac);

                denomenator_wu = MIN(denomenator_gen_jac,MAX(a_sum,b_sum));
                dist_wu = 1.0-numerator_gen_jac/denomenator_wu;

                numerator_sarika = num_sim;
                denomenator_sarika = a_sum+b_sum;
                dist_sarika = 1.0-numerator_sarika/denomenator_sarika;

                rp->normal[idx_out] = dist_jac;
                rp->generalised[idx_out] = dist_gen_jac;
                rp->sarika[idx_out] = dist_sarika;
                rp->wu[idx_out] = dist_wu;
                rp->cosine[idx_out] = result*100;
                idx_out++;
            }
            }
        }

        /* pass out the data */
        //for (i=0;i<idx_out;i++) {
        //    rp->normal[i] = normal[i];
        //    rp->generalised[i] = generalised[i];
        //    rp->sarika[i] = sarika[i];
        //    rp->wu[i] = wu[i];
        //    rp->cosine[i] = cosine[i];
        //}

        //free(normal);
        //free(generalised);
        //free(sarika);
        //free(wu);
        //free(cosine);

        free(data_sum_a);
        free(data_sum_b);
        free(one_data_norm_a);
        free(one_data_norm_b);

        free_dynamic_2d_array_integer(data_jac_a);
        free_dynamic_2d_array_integer(data_jac_b);

        return 0;
    }



    int calc_coeffs_diagnol_triangle(tint **restrict data, tint dim0, tint dim1,
            result_pointers_diagnol *rpd)
    {
        tint i, j, idx_a, idx_b;
        tint *a, *b, *a_jac, *b_jac;
        real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
        real numerator_sarika, denomenator_sarika, dist_sarika;
        real num_sim, numerator_jac, denomenator_jac, numerator_gen_jac, denomenator_gen_jac;
        real a_sum, b_sum, one_an, one_bn, result;

        // calculate some preparation values
        real *restrict data_sum = (real*)malloc(dim0*sizeof(real));
        real *restrict one_data_norm = (real*)malloc(dim0*sizeof(real));
        int **restrict cmbs=NULL;

        int num_cmbs = 0, idx_rpd = 0;
        cmbs = combination_util(dim0,&num_cmbs); 

        //real *restrict normal = (real*)malloc(num_cmbs*sizeof(real));
        //real *restrict generalised = (real*)malloc(num_cmbs*sizeof(real));
        //real *restrict sarika = (real*)malloc(num_cmbs*sizeof(real));
        //real *restrict wu = (real*)malloc(num_cmbs*sizeof(real));
        //real *restrict cosine = (real*)malloc(num_cmbs*sizeof(real));

        tint **restrict data_jac = allocate_dynamic_2d_array_integer(dim0, dim1);

        // preparation values
        //#pragma acc kernels
        for (i=0;i<dim0;i++) {
            data_sum[i] = 0;
            for (j=0;j<dim1;j++) {
                data_jac[i][j]=0;
                if (data[i][j]>0) 
                    data_jac[i][j]=1;
                data_sum[i] += data[i][j];
            }
            //one_data_norm[i] = 1.0/vec_norm(data[i], dim1);
        }

#pragma acc data \
        copy(data[0:dim0][0:dim1],\
                data_sum[0:dim0],\
                data_jac[0:dim0][0:dim1],\
                one_data_norm[0:dim0])
        {

            for (i=0;i< num_cmbs;i++) {
                idx_a = cmbs[i][0], idx_b = cmbs[i][1];

                a = data[idx_a];
                a_sum = data_sum[idx_a];
                a_jac = data_jac[idx_a];
                b = data[idx_b];
                b_sum = data_sum[idx_b];
                b_jac = data_jac[idx_b];

                //vec_add(a, b, summed_array, dim1);
                numerator_jac = sum_minimum_vec(a_jac, b_jac, dim1);

                denomenator_jac = sum_maximum_vec(a_jac, b_jac, dim1);
                //printf("numerator_jac=%f\n",numerator_jac);
                //printf("denomenator_jac=%f\n",denomenator_jac);
                numerator_gen_jac = sum_minimum_vec(a, b, dim1);
                denomenator_gen_jac = sum_maximum_vec(a, b, dim1);

                //printf("numerator_gen_jac=%f\n",numerator_gen_jac);
                //printf("denomenator_gen_jac=%f\n",denomenator_gen_jac);

                num_sim = get_non_zeros_pair(a, b, dim1);
                //printf("num_sim=%f\n",num_sim);

                one_an = one_data_norm[idx_a];
                one_bn = one_data_norm[idx_b];
                real adotb = vec_dot(a,b,dim1);
                //printf("%d %d %.3e\n", idx_a, idx_b, adotb);
                //result = 1.0 - adotb*one_an*one_bn;
                result = adotb*one_an*one_bn;

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
        free_dynamic_2d_array_integer(data_jac);
        free(data_sum);
        free(one_data_norm);


    }

#endif
