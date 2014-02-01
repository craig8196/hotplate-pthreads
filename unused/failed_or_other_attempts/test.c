
#include <stdio.h>
#include <stdlib.h>

#define SIZE 2

int main(int argc, char** argv)
{
    float array[SIZE][SIZE];
    int row, col;
    for(row = 0; row < SIZE; row++)
    {
        for(col = 0; col < SIZE; col++)
        {
            printf("array[%d][%d] at address %#lx or %ld\n", row, col, (long int)&array[row][col], (long int)&array[row][col]);
        }
    }
    return 0;
}
