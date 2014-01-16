/*
 * Each 2D array is organized such that the first index is the row, the second is the column.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


double get_seconds();
void initialize(int size, float plate[size][size]);
void initialize_test_cells(int size, int test[size][size]);
void set_static_cells(int size, float plate[size][size]);
void print_matrix(int size, float plate[size][size]);
int count_cells_by_degrees(int size, float plate[size][size], float temp);
void compute(int size, float* current, float* next, int* test, float error, int* iterations, int* cell_count_gt_50_degrees);


int main(int argc, char* argv[])
{
    int r;
    int repetitions = 1;
    double total_time = 0; // used to eventually calculate an average time
    double fastest_time = FLT_MAX;
    
    for(r = 0; r < repetitions; r++)
    {
        // get start time
        double start_time = get_seconds();
        
        // define counters
        int iteration_count = 0;
        int cell_count_gt_50_degrees = 0;
        
        // define sizes
        const int size = 1024;
        const float error = 0.1f;
        
        // allocate memory to matrices
        float* current_plate = malloc(size*size*sizeof(float));
        float* next_plate = malloc(size*size*sizeof(float));
        int* test = malloc(size*size*sizeof(int));
        
        // check for bad allocation
        if(current_plate == 0 || next_plate == 0)
        {
            printf("Bad allocation.\n");
            return 1;
        }
        
        // initialize the matrices
        initialize(size, (float(*) []) current_plate);
        initialize(size, (float(*) []) next_plate);
        initialize_test_cells(size, (int(*) []) test);
        
        // perform the computation
        compute(size, current_plate, next_plate, test, error, &iteration_count, &cell_count_gt_50_degrees);
        
        // free the matrices
        free(current_plate);
        free(next_plate);
        free(test);
        
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
    
    printf("Average Time: %f\n", total_time/repetitions);
    printf("Fastest Time: %f\n", fastest_time);
    
    return 0;
}

void compute(int size, float* current, float* next, int* test, float error, int* iterations, int* cell_count_gt_50_degrees)
{    
    (*iterations) = 0;
    
    int keep_going = 1;
    
    const int MAX_ITERATIONS = 500;  // a safety while testing
    
#pragma omp parallel shared(keep_going)
    {
        unsigned int it, row, col;
        
        int number_of_threads = omp_get_num_threads();
        int my_thread_num = omp_get_thread_num();
        
        int my_size = size/number_of_threads;
        int start = my_size*my_thread_num;
        int end = start + my_size;
        
        if(start < 1)
            start = 1;
        if(end > size - 1 || my_thread_num == number_of_threads - 1)
            end = size - 1;
        
        // loop to completion
        for(it = 0; it < MAX_ITERATIONS && keep_going; it++)
        {
            // calculate the next iteration
            for(row = start; row < end; row ++)
            {
                float* top = (float*)(current + (row + 1)*size);
                float* curr = (float*)(current + row*size);
                float* bottom = (float*)(current + (row - 1)*size);
                for(col = 1; col < size - 1; col++)
                {
                    (*(next + row*size + col)) = (bottom[col] + top[col] + curr[col - 1] + curr[col + 1] + 4.0f*curr[col])/8.0f;
                }
            }
            
#pragma omp barrier
#pragma omp reduction(||: keep_going)
            {
                keep_going = 0;
                for(row = start; row < end; row++)
                {
                    if(!keep_going)
                    {
                        for(col = 1; col < size - 1; col++)
                        {
                            if(!(*(test + row*size + col)))
                            {
                                float average = ((*(current + (row - 1)*size + col)) +
                                                 (*(current + (row + 1)*size + col)) +
                                                 (*(current + (row)*size + col + 1)) +
                                                 (*(current + (row)*size + col - 1)))/4.0f;
                                           
                                
                                float difference = fabsf((*(current + (row)*size + col)) - average);
                                 
                                if(difference >= error)
                                {
                                    keep_going = 1;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            
#pragma omp barrier
#pragma omp master
            {
                if((*iterations)%50 == 0)
                {
                    printf("Keep going: %d\n", keep_going);
                    print_matrix(size, (float(*)[])current); //*/
                }
                
                // reset certain cells
                set_static_cells(size, (float(*)[])next);
                
                // swap matrices
                float* temp = current;
                current = next;
                next = temp;
                
                // increment the iterations
                if(keep_going)
                {
                    (*iterations)++;
                }
                //printf("Iter: %d\n", (*iterations));
            }
        }
#pragma omp barrier
    }
    
    (*cell_count_gt_50_degrees) = count_cells_by_degrees(size, (float(*) []) next, 50.0f);
    
    return;
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

void initialize_test_cells(int size, int test[size][size])
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




int has_converged(int size, float plate[size][size], float error, int test[size][size])
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
