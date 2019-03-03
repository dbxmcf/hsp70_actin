// Program to print all combination of size r in an array of size n 
#include <stdio.h> 
//#include "ncr.h"
#include "dynamic_2d_array.h"

/* arr[] ---> Input Array 
data[] ---> Temporary array to store current combination 
start & end ---> Staring and Ending indexes in arr[] 
index ---> Current index in data[] 
r ---> Size of a combination to be printed */
void combination_util(int n, int **cmbs) //, int **cmbs, int num_cmb) 
{ 
	int i, j, idx=0;
	// Current combination is ready to be printed, print it 
	for (i=0;i<n;i++)
	{ 
        for (j=i+1; j<n; j++) {
			cmbs[idx][0]=i;
            cmbs[idx][1]=j;
			idx++;
		}
	} 

} 

// Driver program to test above functions 
int main() 
{ 
	int n=5,r=2; 
	int **cmbs=NULL;
    int num_cmbs = n*(n-1)/2;
    printf("num_cmbs=%d\n",num_cmbs);
    cmbs = allocate_dynamic_2d_array(num_cmbs,2);
    printf("n=%d\n",n);
	combination_util(n,cmbs); 
    print_matrix(cmbs, num_cmbs, r, "%3d");
    free_dynamic_2d_array(cmbs);
} 
