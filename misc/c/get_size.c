#include <stdio.h>
#include <stdbool.h>
#include <limits.h>

int main() {
    printf("Storage size for int : %d \n", sizeof(int));
    printf("Storage size for short : %d \n", sizeof(short));
    printf("Storage size for char : %d \n", sizeof(char));
    printf("Storage size for bool : %d \n", sizeof(bool));
    return 0;
}
