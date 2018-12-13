#include <stdio.h>

#define N 2
#define D 5
#define C 3

typedef float real;
typedef uint8 tint;
// https://www.geeksforgeeks.org/compute-the-minimum-or-maximum-max-of-two-integers-without-branching/
// https://stackoverflow.com/questions/24529504/find-out-max-min-of-two-number-without-using-if-else
// http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
tint get_minimum_vec(tint *a, tint, *b, tint *c, int vec_dim)
{
    for (int i=0;i<vec_dim;i++)
        c[i] = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // min(x, y)
}

tint get_maximum_vec(tint *a, tint, *b, tint *c, int vec_dim)
{
    for (int i=0;i<vec_dim;i++)
        c[i] = b[i] ^ ((a[i] ^ b[i]) & -(a[i] < b[i])); // max(x, y)
}

tint get_non_zeros(tint *a, tint, *b, tint *c, int vec_dim)
{
    for (int i=0;i<vec_dim;i++)
        
}

//real get_sum_vec(real a, real b, real c, int vec_dim)
//{
//    real c[];
//    for (int i=0;i<vec_dim;i++)
//        if (a[i] > b[i]) c[i] = a[i]+b[i];
//}
//
//real vec_add(real a, real b, real c, int vec_dim)
//{
//    for(i=0; i<vec_dim; i++) 
//        c[i] = a[i]+b[i];
//}

int main(void)
{
    
    tint data[N][D] = {{1,2,3,4,5},{6,7,8,9,10},{7,8,9,10,11}};
    tint cmb[C][2]={{0,1},{0,2},{1,2}}
    float data_sum[N];
    int data_jac[N][D]={0};

    for (int i=0;i<N;i++)
        data_sum[i] = 0;
        for (int j=0;j<D;j++)
            if (data[i][j]>0) data_jac[i][j]=1;
            data_sum[i] += data[i][j]

    int idx_a, idx_b;
    int non_zeros[D]={0};
    
    for (int i=0;i<N;i++) {
        idx_a = cmb[i][0], idx_b = cmb[i][1];

        a = data[idx_a]
        a_sum = data_sum[idx_a]
        a_jac = data_jac[idx_a]
        b = data[idx_b]
        b_sum = data_sum[idx_b]
        b_jac = data_jac[idx_b]

        for (int j=0;j<D;j++)
            
        non_zeros = (a >0) & (b > 0)
        summed_array = a + b

        numerator_jac = np.sum(np.minimum(a_jac,b_jac))
        denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
        numerator_gen_jac =np.sum(np.minimum(a,b))
        denomenator_gen_jac =np.sum(np.maximum(a,b))
        num_sim = np.sum(summed_array[non_zeros])
        result = 1 - spatial.distance.cosine(a, b)
    }
}
