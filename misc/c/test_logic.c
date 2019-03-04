#include <stdio.h>

int main() {

    short a = 5;
    short b = 3;
    short c = 0;

    c= a && b;
    printf("a && b = %d\n", c);
    c= a || b;
    printf("a || b = %d\n", c);

    a = 5, b = 0;
    c= a && b;
    printf("a && b = %d\n", c);
    c= a || b;
    printf("a || b = %d\n", c);

}
