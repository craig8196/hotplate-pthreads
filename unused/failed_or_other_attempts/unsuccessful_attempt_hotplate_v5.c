/*
 * Each 2D array is organized such that the first index is the row, the second is the column.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

double get_seconds();
void initialize(int size, float plate[size][size]);
void initialize_test_cells(int size, char test[size][size]);
void set_static_cells(int size, float plate[size][size]);
void swap(float** current, float** next);
void print_matrix(int size, float plate[size][size]);
int count_cells_by_degrees(int size, float plate[size][size], float temp);
void hotplate(int size, float error, int* iterations, int* cell_count_gt_50_degrees);
inline void compute(int size, float current[size][size], float next[size][size]);
inline int has_converged(int size, float plate[size][size], float error, char test[size][size]);

int main(int argc, char* argv[])
{
    int r;
    int repetitions = 10;
    double total_time = 0;
    double fastest_time = FLT_MAX;
    
    for(r = 0; r < repetitions; r++)
    {
        // get start time
        double start_time = get_seconds();

        
        // define sizes
        const int size = 1024;
        const float error = 0.1f;
        // define counters
        int iteration_count = 0;
        int cell_count_gt_50_degrees = 0;
        // run hotplate
        hotplate(size, error, &iteration_count, &cell_count_gt_50_degrees);
        
        
        // get stop time
        double end_time = get_seconds();
        double time_interval = end_time - start_time;
        
        // report convergence and time
        printf("Iterations: %d\n", iteration_count);
        printf("Cells with >= 50.0 degrees: %d\n", cell_count_gt_50_degrees);
        printf("%d %f\n", omp_get_max_threads(), time_interval); // number_of_threads time_to_execute
        fflush(stdout);
        
        total_time += time_interval;
        if(time_interval < fastest_time)
        {
            fastest_time = time_interval;
        }
    }
    
    // report average and fastest times
    printf("Average Time: %f\n", total_time/repetitions);
    printf("Fastest Time: %f\n", fastest_time);
    fflush(stdout);
    
    return 0;
}

void hotplate(int size, float error, int* iterations, int* cell_count_gt_50_degrees)
{
    float** matrices;
    
    // allocate memory to matrices
    char* test = malloc(size*size*sizeof(char));
    
    // initialize the matrices
    initialize_test_cells(size, (char(*) []) test);
    
    (*iterations) = 0;
    
    int keep_going = 1;
    int* keep_goings;
    
    const int MAX_ITERATIONS = 500;  // a safety while testing
    
    #pragma omp parallel shared(keep_going, keep_goings, matrices, cell_count_gt_50_degrees, test, iterations)
    {
        unsigned int it, nt, row, col;
        int num_threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();
#pragma omp single
        {
            matrices = malloc(sizeof(float*)*(num_threads+1));
            for(nt = 0; nt < num_threads + 1; nt++)
            {
                matrices[nt] = malloc(size*size*sizeof(float));
                initialize(size, (float(*) []) matrices[nt]);
            }
            keep_goings = malloc(sizeof(int)*(num_threads));
        }
        
#pragma omp barrier
        
        // loop to completion
        for(it = 0; it < MAX_ITERATIONS && keep_going; it+=num_threads)
        {
                
            // calculate the next several iterations
            for(nt = 0; nt < num_threads; nt++)
            {
                compute(size, (float(*) [])matrices[nt], (float(*) [])matrices[nt+1]);
                set_static_cells(size, (float(*)[])matrices[thread_num+1]);
#pragma omp barrier
            }
            
            
            
            keep_goings[thread_num] = has_converged(size, (float(*)[])matrices[thread_num], error, (char(*)[]) test);

#pragma omp barrier
#pragma omp master
            {
                // check if we should keep going
                keep_going = 1;
                for(nt = 0; nt < num_threads; nt++)
                {
                    keep_going = keep_going && !keep_goings[nt];
                    if(!keep_goings[nt])
                    {
                        (*iterations)++;
                    }
                    else
                    {
                        (*cell_count_gt_50_degrees) = count_cells_by_degrees(size, (float(*) []) matrices[nt], 50.0f);
                        break;
                    }
                }
                
                // swaps
                float* temp = matrices[0];
                matrices[0] = matrices[num_threads];
                matrices[num_threads] = temp;
            }
#pragma omp barrier
        }
        
        
#pragma omp single
        {
            for(nt = 0; nt < num_threads + 1; nt++)
            {
                free(matrices[nt]);
            }
            free(keep_goings);
            free(matrices);
        }

    }
    
    free(test);
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

void initialize_test_cells(int size, char test[size][size])
{
    int row, col;
    
    for(row = 0; row < size; row++)
    {
        for(col = 0; col < size; col++)
        {
            test[row][col] = 0;
        }
    }
    
    if(size > 500)
    {
        for(col = 0; col < 331; col++)
        {
            test[400][col] = 1;
        }
    
        test[200][500] = 1;
    }
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

inline int has_converged(int size, float plate[size][size], float error, char test[size][size])
{
    int row, col;
    int converged = 1;

    for(row = 1; row < size - 1; row++)
    {
        for(col = 1; col < size - 1; col++)
        {
            if(!test[row][col])
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
        if(!converged)
        {
            break;
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

void swap(float** current, float** next)
{
    float* temp = *current;
    *current = *next;
    *next = temp;
}

inline void compute(int size, float current[size][size], float next[size][size])
{
    int row, col;
    
#pragma omp for schedule(dynamic, 64) private(row, col) nowait
    for(row = 1; row < size - 1; row++)
    {
        float* top = current[row+1];
        float* curr = current[row];
        float* bottom = current[row-1];
        for(col = 1; col < size - 1; col++)
        {
            next[row][col] = (bottom[col] + top[col] + curr[col - 1] + curr[col + 1] + 4.0f*curr[col])/8.0f;
        }
    }
}
