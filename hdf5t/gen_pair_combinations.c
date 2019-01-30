// Program to print all combination of size r in an array of size n 
#include <stdio.h> 
#include "ncr.h"
#include "dynamic_2d_array.h"


// The main function that prints all combinations of size r 
// in arr[] of size n. This function mainly uses combination_util() 
//void printCombination(int arr[], int n, int r) 
//{ 
//	// A temporary array to store all combination one by one 
//	int i,data[r],idx;
//	int **cmbs=NULL;
//	int num_cmb=0, idx_cmb=0;
//	num_cmb = ncr(5,3);
//	printf("num_cmb=%d\n",num_cmb);
//	// Print all combination using temprary array 'data[]' 
//	cmbs = allocate_dynamic_2d_array(num_cmb,2);
//	combination_util(arr, data, 0, n-1, 0, r, cmbs, num_cmb, &idx_cmb);
//	free_dynamic_2d_array(cmbs);
//} 

/* arr[] ---> Input Array 
data[] ---> Temporary array to store current combination 
start & end ---> Staring and Ending indexes in arr[] 
index ---> Current index in data[] 
r ---> Size of a combination to be printed */
void combination_util(int n, int **cmbs) //, int **cmbs, int num_cmb) 
{ 
	int i, j, idx;
	// Current combination is ready to be printed, print it 
	for (i=0;i<n;i++)
	{ 
        for (j=i+1; j<n; j++) {
			//printf("%d %d\n", i, j); 
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
