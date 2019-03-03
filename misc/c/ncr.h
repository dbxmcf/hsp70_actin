#ifndef _NCR_
#define _NCR_
/*
 * C program to Calculate the value of nCr
 */
#include <stdio.h>
 
int fact(int z);
 
int ncr(int n, int r)
{
    //printf("\n Enter the value for N and R \n");
    //scanf("%d%d", &n, &r);
    int val = fact(n) / (fact(r) * fact(n - r));
    //printf("\n The value of ncr is: %d", ncr);
    return val;
}
 
int fact(int z)
{
    int f = 1, i;
    if (z == 0)
    {
        return(f);
    }
    else
    {
        for (i = 1; i <= z; i++)
	{
            f = f * i;
	}
    }
    return(f);
}
#endif
