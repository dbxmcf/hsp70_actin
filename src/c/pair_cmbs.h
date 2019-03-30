// Program to print all combination of size r in an array of size n 
#include <stdio.h> 
//#include "ncr.h"
#include "dynamic_2d_array.h"

/* arr[] ---> Input Array 
data[] ---> Temporary array to store current combination 
start & end ---> Staring and Ending indexes in arr[] 
index ---> Current index in data[] 
r ---> Size of a combination to be printed */
integer **combination_util(integer n, integer *num_cmbs) //, int **cmbs, int num_cmb) 
{ 
	integer i, j, idx=0;
	*num_cmbs = n*(n-1)/2;
	integer **cmbs = allocate_dynamic_2d_array_integer(*num_cmbs,2);
	// Current combination is ready to be printed, print it 
	for (i=0;i<n;i++)
	{ 
        for (j=i+1; j<n; j++) {
			cmbs[idx][0]=i;
            cmbs[idx][1]=j;
			idx++;
		}
	} 
	return cmbs;
} 

void distribute_parts_start_size(integer part_dim, integer num_parts, integer **ptr_part_start, integer **ptr_part_size)
{
	integer i;
	integer *part_start = (integer*)malloc(num_parts*sizeof(integer));
	integer *part_size  = (integer*)malloc(num_parts*sizeof(integer));
	integer part_avg_lines = part_dim / num_parts;
    integer part_rmd_lines = part_dim % num_parts;

    //dvc_blk_part_a_start[0] = 0;
    for (i=0;i<part_rmd_lines;i++) {
		part_size[i] = part_avg_lines + 1;
		//printf("part_size[%ld]-%ld\n",i,part_size[i]);
    }
    for (i=part_rmd_lines;i<num_parts;i++) {
		part_size[i] = part_avg_lines;
		//printf("part_size[%ld]-%ld\n",i,part_size[i]);
    }

	part_start[0] = 0;
    for (i=1;i<num_parts;i++) {
		part_start[i] = part_start[i-1] + part_size[i-1];
		//printf("part_start[%ld]-%ld\n",i,part_start[i]);
    }
	*ptr_part_start = part_start;
	*ptr_part_size = part_size;
    //printf("here---%d\n",mpi_rank);
    //printf("here---\n");
}

// Driver program to test above functions 
//int main() 
//{ 
//	int n=5,r=2; 
//	int **cmbs=NULL;
//    int num_cmbs = n*(n-1)/2;
//    printf("num_cmbs=%d\n",num_cmbs);
//    cmbs = allocate_dynamic_2d_array(num_cmbs,2);
//    printf("n=%d\n",n);
//	combination_util(n,cmbs); 
//    print_matrix(cmbs, num_cmbs, r, "%3d");
//    free_dynamic_2d_array(cmbs);
//} 
