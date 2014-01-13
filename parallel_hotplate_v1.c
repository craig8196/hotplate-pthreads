/*
 * Each 2D array is organized such that the first index is the row, the second is the column.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 1024
#define ERROR 0.1

double get_seconds();
void initialize(int size, float plate[size][size]);
void set_static_cells(int size, float plate[size][size]);
void swap(float** current, float** next);
void next_iteration(int size, float current[size][size], float next[size][size]);
int has_converged(int size, float current[size][size], float error);
void print_matrix(int size, float plate[size][size]);
int count_cells_by_degrees(int size, float plate[size][size], float temp);

int main(int argc, char *argv[])
{
    // get start time
    double start_time = get_seconds();
    
    int iteration_count = 0;
    
    // allocate memory to matrices
    float* current_plate = malloc(SIZE*SIZE*sizeof(float));
    float* next_plate = malloc(SIZE*SIZE*sizeof(float));
    
    // check for bad allocation
    if(current_plate == 0 || next_plate == 0)
    {
        printf("Bad allocation.\n");
        return 1;
    }
    
    // initialize the matrices
    initialize(SIZE, (float(*) [SIZE]) current_plate);
    initialize(SIZE, (float(*) [SIZE]) next_plate);
    
    // iterate to convergence
    while(!has_converged(SIZE, (float(*) [SIZE]) current_plate, ERROR))
    {
        next_iteration(SIZE, (float(*) [SIZE]) current_plate, (float(*) [SIZE]) next_plate);
        swap(&current_plate, &next_plate);
        iteration_count++;
    }
    
    int cell_count_gt_50_degrees = count_cells_by_degrees(SIZE, (float(*) [SIZE]) current_plate, 50.0f);
    
    // free the matrices
    free(current_plate);
    free(next_plate);
    
    // get stop time
    double end_time = get_seconds();
    double total_time = end_time - start_time;
    
    // report convergence and time
    printf("Iterations: %d\n", iteration_count);
    printf("Time: %f\n", total_time);
    printf("Cells with >= 50.0 degrees: %d\n", cell_count_gt_50_degrees);
    fflush(stdout);
    
    return 0;
}

double get_seconds()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

void initialize(int size, float plate[size][size])
{
    int row, col;
    
    // initialize left and right sides
    for(row = 0; row < size; row++)
    {
        plate[row][0] = 0.0;
        plate[row][size - 1] = 0.0;
    }
    
    // initialize top and bottom
    for(col = 0; col < size; col++)
    {
        plate[0][col] = 0.0;
    }
    
    for(col = 0; col < size; col++)
    {
        plate[size - 1][col] = 100.0;
    }
    
    // initialize center area
    for(row = 1; row < size - 1; row++)
    {
        for(col = 1; col < size - 1; col++)
        {
            plate[row][col] = 50.0;
        }
    }
    
    // initialize special cells
    set_static_cells(size, plate);
}

void set_static_cells(int size, float plate[size][size])
{
    int col;
    
    if(size > 500)
    {
        for(col = 0; col < 331; col++)
        {
            plate[400][col] = 100.0;
        }
    
        plate[200][500] = 100.0;
    }
}


void swap(float** current, float** next)
{
    float* temp = *current;
    *current = *next;
    *next = temp;
}


void next_iteration(int size, float current[size][size], float next[size][size])
{
    int row, col;
    
    #pragma omp parallel for private(col)
    for(row = 1; row < size - 1; row++)
    {
        for(col = 1; col < size - 1; col++)
        {
            next[row][col] = (current[row - 1][col] + current[row + 1][col] + 
                              current[row][col - 1] + current[row][col + 1] + 
                              4.0*current[row][col])/8.0f; 
        }
    }
    
    set_static_cells(size, next);
}

int has_converged(int size, float plate[size][size], float error)
{
    int row, col;
    int converged = 1;
    
    #pragma omp parallel for private(col)
    for(row = 1; row < size - 1; row++)
    {
        if(converged)
            for(col = 1; col < size - 1; col++)
            {
                if(!((row == 400 && col > 0 && col <= 330) || (row == 200 && col == 500)))
                {
                    float average = (plate[row - 1][col] + plate[row + 1][col] + 
                                     plate[row][col - 1] + plate[row][col + 1])/4.0f;
                               
                    
                    float difference = fabsf(plate[row][col] - average);
                    
                    /*printf("Avg: %f", average);
                    printf("Val: %f", plate[row][col]);
                    printf("Dif: %f", difference);*/
                        
                    if(difference >= error)
                    {
                        converged = 0;
                        break;
                    }
                }
            }
    }
    
    return converged;
}

void print_matrix(int size, float plate[size][size])
{
    int row, col;
    for(row = 0; row < size && row < 10; row++)
    {
        for(col = 0; col < size && col < 10; col++)
        {
            printf("V%f ", plate[row][col]);
        }
        printf("\n");
    }
    printf("\n");
}

int count_cells_by_degrees(int size, float plate[size][size], float temp)
{
    int count = 0;
    int row, col;
    
    for(row = 1; row < size - 1; row++)
    {
        for(col = 1; col < size - 1; col++)
        {
            if(plate[row][col] >= temp)
            {
                count++;
            }
        }
    }
    
    return count;
}
