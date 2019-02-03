#ifndef _CALC_COEFFS_H_
#define _CALC_COEFFS_H_
#include <stdio.h>
#include <math.h>

#define N 3
#define D 5
#define C 3

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


// https://www.geeksforgeeks.org/compute-the-minimum-or-maximum-max-of-two-integers-without-branching/
// https://stackoverflow.com/questions/24529504/find-out-max-min-of-two-number-without-using-if-else
// http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
real sum_minimum_vec(tint *a, tint *b, tint vec_dim)
{
  tint i, c;  
  real sum=0;
  for (i=0;i<vec_dim;i++) {
        // c = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // min(x, y)
        c = ((a[i])<(b[i]))?(a[i]):(b[i]);
        //printf("%d, %d, %d\n", a[i],b[i],c);
        sum += c;
  }
  return sum;
}

real sum_maximum_vec(tint *a, tint *b, tint vec_dim)
{
  tint i, c;  
  real sum=0;  
  for (i=0;i<vec_dim;i++) {
        //c = a[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // max(x, y)
        c = ((a[i])>(b[i]))?(a[i]):(b[i]);
        sum += c;
  }
  return sum;
}

real get_non_zeros_pair(tint *a, tint *b, tint vec_dim)
{
    real sum = 0;
    tint i;
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

real vec_norm(tint *a, tint vec_dim)
{
    tint i;
    real norm=0;
    for(i=0;i<vec_dim; i++) {
        norm += a[i]*a[i];
    }
    return sqrt(norm);
}

real vec_dot(tint *a, tint *b, tint vec_dim)
{
    tint i;
    real sum=0;
    for(i=0;i<vec_dim; i++) {
        sum += a[i]*b[i];
    }
    return sum;
}

// https://stackoverflow.com/questions/20013693/read-csv-file-to-a-2d-array-on-c
void read_h5()
{
    
}

int calc_coeffs_block(tint **data_part_a, tint part_a_dim0, tint part_a_dim1,
                      tint **data_part_b, tint part_b_dim0, tint part_b_dim1,
                      result_pointers *rp)
{
    //int i, j, idx_a, idx_b;
    tint i, j, idx_a, idx_b;
    tint *a, *b, *a_jac, *b_jac;
    real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
    real numerator_sarika, denomenator_sarika, dist_sarika;
    real num_sim, numerator_jac, denomenator_jac, numerator_gen_jac, denomenator_gen_jac;

    // perform full permutation of all pairs of part a and part b
    if (part_a_dim1 != part_b_dim1) {
        printf("dim1s do not match!");
        exit(1);
    }

    tint dim1 = part_a_dim1;
    // calculate some preparation values
    real *data_sum_a = (real*)malloc(part_a_dim0*sizeof(real));
    real *data_sum_b = (real*)malloc(part_b_dim0*sizeof(real));
    real *one_data_norm_a = (real*)malloc(part_a_dim0*sizeof(real));
    real *one_data_norm_b = (real*)malloc(part_b_dim0*sizeof(real));
    real a_sum, b_sum, one_an, one_bn, result;

    tint **data_jac_a = allocate_dynamic_2d_array_integer(part_a_dim0, dim1);
    tint **data_jac_b = allocate_dynamic_2d_array_integer(part_b_dim0, dim1);

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
    for (i=0;i<part_a_dim0;i++) {
        data_sum_a[i] = 0;
        for (j=0;j<part_a_dim1;j++) {
            if (data_part_a[i][j]>0) 
                data_jac_a[i][j]=1;
            data_sum_a[i] += data_part_a[i][j];
        }
        one_data_norm_a[i] = 1.0/vec_norm(data_part_a[i], dim1);
    }

    // part_b values
    for (i=0;i<part_b_dim0;i++) {
        data_sum_b[i] = 0;
        for (j=0;j<part_b_dim1;j++) {
            if (data_part_b[i][j]>0) 
                data_jac_b[i][j]=1;
            data_sum_b[i] += data_part_b[i][j];
        }
        one_data_norm_b[i] = 1.0/vec_norm(data_part_b[i], dim1);
    }

    // the large loop that calculates the matrix
    for (idx_a=0;idx_a<part_a_dim0;idx_a++) {
        for (idx_b=0;idx_b<part_b_dim0;idx_b++){
            // 
            a = data_part_a[idx_a];
            a_sum = data_sum_a[idx_a];
            a_jac = data_jac_a[idx_a];
            b = data_part_b[idx_b];
            b_sum = data_sum_b[idx_b];
            b_jac = data_jac_b[idx_b];

            //vec_add(a, b, summed_array, D);
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
            result = 1.0 - vec_dot(a,b,dim1)*one_an*one_bn;
            
            dist_gen_jac = 1.0-numerator_gen_jac/denomenator_gen_jac;
            dist_jac = 1.0-numerator_jac/denomenator_jac;
            
            denomenator_wu = MIN(denomenator_gen_jac,MAX(a_sum,b_sum));
            dist_wu = 1.0-numerator_gen_jac/denomenator_wu;
            
            numerator_sarika = num_sim;
            denomenator_sarika = a_sum+b_sum;
            dist_sarika = 1.0-numerator_sarika/denomenator_sarika;

            rp->normal[idx_a][idx_b] = dist_jac;
            rp->generalised[idx_a][idx_b] = dist_gen_jac;
            rp->sarika[idx_a][idx_b] = dist_sarika;
            rp->wu[idx_a][idx_b] = dist_wu;
            rp->cosine[idx_a][idx_b] = result*100;
        }
    }

    free(data_sum_a);
    free(data_sum_b);
    free(one_data_norm_a);
    free(one_data_norm_b);

    free_dynamic_2d_array_integer(data_jac_a);
    free_dynamic_2d_array_integer(data_jac_b);

    //print_matrix_real(normal, part_a_dim0, part_b_dim0, "%5.3f ");
    //free_dynamic_2d_array_integer(cmb_ab);

    ////tint data[N][D] = {{1,2,3,4,5},{6,7,8,9,10},{7,8,9,10,11}};
    ////tint data[N][D] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    //tint data[N][D] = {{1,0,3,0,5},{6,7,0,9,10},{11,0,13,14,0}};
    //tint data_jac[N][D]={0};
    //tint cmb[C][2]={{0,1},{0,2},{1,2}};
    //real data_sum[N], a_sum, b_sum, one_data_norm[N], one_an, one_bn, result;
    //real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
    //real numerator_sarika, denomenator_sarika, dist_sarika;
//
    //tint i, j, idx_a, idx_b;
    //tint *a, *b, *a_jac, *b_jac;
    ////tint a[D], b[D], a_jac[D], b_jac[D];
    //tint non_zeros[D]={0};
    //tint summed_array[D]={0};
    //real num_sim, numerator_jac, denomenator_jac, numerator_gen_jac, denomenator_gen_jac;
    //    
    //for (i=0;i<N;i++) {
    //    data_sum[i] = 0;
    //    for (j=0;j<D;j++) {
    //        if (data[i][j]>0) 
    //            data_jac[i][j]=1;
    //        data_sum[i] += data[i][j];
    //    }
    //    one_data_norm[i] = 1.0/vec_norm(data[i], D);
    //}
//
    //for (i=0;i< C;i++) {
    //    idx_a = cmb[i][0], idx_b = cmb[i][1];
//
    //    a = data[idx_a];
    //    a_sum = data_sum[idx_a];
    //    a_jac = data_jac[idx_a];
    //    b = data[idx_b];
    //    b_sum = data_sum[idx_b];
    //    b_jac = data_jac[idx_b];
//
    //    //vec_add(a, b, summed_array, D);
    //    numerator_jac = sum_minimum_vec(a_jac, b_jac, D);
//
    //    denomenator_jac = sum_maximum_vec(a_jac, b_jac, D);
    //    //printf("numerator_jac=%f\n",numerator_jac);
    //    //printf("denomenator_jac=%f\n",denomenator_jac);
    //    numerator_gen_jac = sum_minimum_vec(a, b, D);
    //    denomenator_gen_jac = sum_maximum_vec(a, b, D);
    //    
    //    //printf("numerator_gen_jac=%f\n",numerator_gen_jac);
    //    //printf("denomenator_gen_jac=%f\n",denomenator_gen_jac);
    //    
    //    num_sim = get_non_zeros_pair(a, b, D);
    //    //printf("num_sim=%f\n",num_sim);
    //    
    //    one_an = one_data_norm[idx_a];
    //    one_bn = one_data_norm[idx_b];
    //    result = 1.0 - vec_dot(a,b,D)*one_an*one_bn;
    //    
    //    dist_gen_jac = 1.0-numerator_gen_jac/denomenator_gen_jac;
    //    dist_jac = 1.0-numerator_jac/denomenator_jac;
    //    
    //    denomenator_wu = MIN(denomenator_gen_jac,MAX(a_sum,b_sum));
    //    dist_wu = 1.0-numerator_gen_jac/denomenator_wu;
    //    
    //    numerator_sarika = num_sim;
    //    denomenator_sarika = a_sum+b_sum;
    //    dist_sarika = 1.0-numerator_sarika/denomenator_sarika;
    //    
    //    //printf("normal = %f\n", dist_jac);
    //    //printf("generalised =%f\n", dist_gen_jac);
    //    //printf("sarika =%f\n", dist_sarika);
    //    //printf("wu =%f\n", dist_wu);
    //    //printf("result = %f\n", result);
    //    
    //    //non_zeros = (a >0) & (b > 0)
    //    //summed_array = a + b
    //    //
    //    //numerator_jac = np.sum(np.minimum(a_jac,b_jac))
    //    //denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
    //    //numerator_gen_jac =np.sum(np.minimum(a,b))
    //    //denomenator_gen_jac =np.sum(np.maximum(a,b))
    //    //num_sim = np.sum(summed_array[non_zeros])
    //    //result = 1 - spatial.distance.cosine(a, b)
    //    
    //    //dist_gen_jac=1.0-(float(numerator_gen_jac)/float(denomenator_gen_jac))                    
    //    //dist_jac=1.0-(float(numerator_jac)/float(denomenator_jac))
    //    //
    //    //denomenator_wu = min(denomenator_gen_jac,max(a_sum,b_sum) )
    //    //dist_wu = 1.0-(float(numerator_gen_jac)/float(denomenator_wu))
    //    //
    //    //numerator_sarika = num_sim
    //    //denomenator_sarika = a_sum+b_sum
    //    //dist_sarika = 1.0-(float(numerator_sarika)/float(denomenator_sarika))
    
}

#endif