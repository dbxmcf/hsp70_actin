#include <stdio.h>
#include <math.h>

#define N 17
#define D 12
#define C 3

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


typedef float real;
typedef unsigned short tint;

// https://www.geeksforgeeks.org/compute-the-minimum-or-maximum-max-of-two-integers-without-branching/
// https://stackoverflow.com/questions/24529504/find-out-max-min-of-two-number-without-using-if-else
// http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
real sum_minimum_vec(tint *a, tint *b, tint vec_dim)
{
  tint i, c;  
  real sum=0;
  for (i=0;i<vec_dim;i++) {
        c = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // min(x, y)
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
        c = a[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // max(x, y)
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

int main(void)
{
    
    //tint data[N][D] = {{1,2,3,4,5},{6,7,8,9,10},{7,8,9,10,11}};
    //tint data[N][D] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    //tint data[N][D] = {{1,0,3,0,5},{6,7,0,9,10},{11,0,13,14,0}};
    tint data[N][D] = {{  0,    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11},
                       {100,  101,  102,  103,  104,  105,  106,  107,  108,  109,  110,  111},
                       {200,  201,  202,  203,  204,  205,  206,  207,  208,  209,  210,  211},
                       {300,  301,  302,  303,  304,  305,  306,  307,  308,  309,  310,  311},
                       {400,  401,  402,  403,  404,  405,  406,  407,  408,  409,  410,  411},
                       {500,  501,  502,  503,  504,  505,  506,  507,  508,  509,  510,  511},
                       {600,  601,  602,  603,  604,  605,  606,  607,  608,  609,  610,  611},
                       {700,  701,  702,  703,  704,  705,  706,  707,  708,  709,  710,  711},
                       {800,  801,  802,  803,  804,  805,  806,  807,  808,  809,  810,  811},
                       {900,  901,  902,  903,  904,  905,  906,  907,  908,  909,  910,  911},
                       {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011},
                       {1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111},
                       {1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211},
                       {1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308, 1309, 1310, 1311},
                       {1400, 1401, 1402, 1403, 1404, 1405, 1406, 1407, 1408, 1409, 1410, 1411},
                       {1500, 1501, 1502, 1503, 1504, 1505, 1506, 1507, 1508, 1509, 1510, 1511},
                       {1600, 1601, 1602, 1603, 1604, 1605, 1606, 1607, 1608, 1609, 1610, 1611}};
    tint data_jac[N][D]={0};
    tint cmb[C][2]={{0,1},{0,2},{1,2}};
    real data_sum[N], a_sum, b_sum, one_data_norm[N], one_an, one_bn, result;
    real dist_gen_jac, dist_jac, denomenator_wu, dist_wu; 
    real numerator_sarika, denomenator_sarika, dist_sarika;

    tint i, j, idx_a, idx_b;
    tint *a, *b, *a_jac, *b_jac;
    //tint a[D], b[D], a_jac[D], b_jac[D];
    tint non_zeros[D]={0};
    tint summed_array[D]={0};
    real num_sim, numerator_jac, denomenator_jac, numerator_gen_jac, denomenator_gen_jac;
        
    for (i=0;i<N;i++) {
        data_sum[i] = 0;
        for (j=0;j<D;j++) {
            if (data[i][j]>0) 
                data_jac[i][j]=1;
            data_sum[i] += data[i][j];
        }
        one_data_norm[i] = 1.0/vec_norm(data[i], D);
    }

    for (i=0;i< C;i++) {
        idx_a = cmb[i][0], idx_b = cmb[i][1];

        a = data[idx_a];
        a_sum = data_sum[idx_a];
        a_jac = data_jac[idx_a];
        b = data[idx_b];
        b_sum = data_sum[idx_b];
        b_jac = data_jac[idx_b];

        //vec_add(a, b, summed_array, D);
        numerator_jac = sum_minimum_vec(a_jac, b_jac, D);

        denomenator_jac = sum_maximum_vec(a_jac, b_jac, D);
        //printf("numerator_jac=%f\n",numerator_jac);
        //printf("denomenator_jac=%f\n",denomenator_jac);
        numerator_gen_jac = sum_minimum_vec(a, b, D);
        denomenator_gen_jac = sum_maximum_vec(a, b, D);
        
        //printf("numerator_gen_jac=%f\n",numerator_gen_jac);
        //printf("denomenator_gen_jac=%f\n",denomenator_gen_jac);
        
        num_sim = get_non_zeros_pair(a, b, D);
        //printf("num_sim=%f\n",num_sim);
        
        one_an = one_data_norm[idx_a];
        one_bn = one_data_norm[idx_b];
        result = 1.0 - vec_dot(a,b,D)*one_an*one_bn;
        
        dist_gen_jac = 1.0-numerator_gen_jac/denomenator_gen_jac;
        dist_jac = 1.0-numerator_jac/denomenator_jac;
        
        denomenator_wu = MIN(denomenator_gen_jac,MAX(a_sum,b_sum));
        dist_wu = 1.0-numerator_gen_jac/denomenator_wu;
        
        numerator_sarika = num_sim;
        denomenator_sarika = a_sum+b_sum;
        dist_sarika = 1.0-numerator_sarika/denomenator_sarika;
        
        printf("normal = %f\n", dist_jac);
        printf("generalised =%f\n", dist_gen_jac);
        printf("sarika =%f\n", dist_sarika);
        printf("wu =%f\n", dist_wu);
        printf("result = %f\n", result*100);
        
        //non_zeros = (a >0) & (b > 0)
        //summed_array = a + b
        //
        //numerator_jac = np.sum(np.minimum(a_jac,b_jac))
        //denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
        //numerator_gen_jac =np.sum(np.minimum(a,b))
        //denomenator_gen_jac =np.sum(np.maximum(a,b))
        //num_sim = np.sum(summed_array[non_zeros])
        //result = 1 - spatial.distance.cosine(a, b)
        
        //dist_gen_jac=1.0-(float(numerator_gen_jac)/float(denomenator_gen_jac))                    
        //dist_jac=1.0-(float(numerator_jac)/float(denomenator_jac))
        //
        //denomenator_wu = min(denomenator_gen_jac,max(a_sum,b_sum) )
        //dist_wu = 1.0-(float(numerator_gen_jac)/float(denomenator_wu))
        //
        //numerator_sarika = num_sim
        //denomenator_sarika = a_sum+b_sum
        //dist_sarika = 1.0-(float(numerator_sarika)/float(denomenator_sarika))
    }
}
