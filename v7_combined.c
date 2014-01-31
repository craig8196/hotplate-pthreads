/*
 * Each 2D array is organized such that the first index is the row, the second is the column.
 * 
 * Version 6: Dynamic scheduler.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <float.h>

// structs
typedef struct thread_info_t
{
    int my_thread_num;
    int num_of_threads;
    
    int size;
    float** current_hotplate;
    float** next_hotplate;
    int chunk_size;
    int* chunk_index;
    int chunk_count;
    pthread_mutex_t* chunk_lock;
    
    float error;
    int** test_plate;
    int* keep_going;
    int* row_index;
    int* col_index;
    
    int* iterations;
    int* cell_count;
    double my_compute_time;
    double my_converged_time;
    int my_converged_count;
    double my_wait_time;
    
    int* count;
    pthread_mutex_t* count_lock;
    volatile int my_wait;
    
    struct thread_info_t* all_threads_info;
    
    char cache_buffer[64];
} ThreadInfo;



// function signatures
double get_seconds();
void initialize(int size, float plate[size][size]);
void initialize_test_cells(int size, int test[size][size]);
void set_static_cells(int size, float plate[size][size]);
void swap(float** current, float** next);
void compute(int next_chunk_index, int chunk_size, int size, float current[size][size], float next[size][size]);
int quick_convergence_test(ThreadInfo* info, int size, float plate[size][size], float error, int test[size][size]);
void* run_thread(void* thread_info);
void run_hotplate(int num_threads, int size_of_plate, float error, int* iterations, int* cell_count_gt_50_degrees);
void barrier(ThreadInfo* info);
void compute_static(ThreadInfo* info, int size, float current[size][size], float next[size][size]);


// performs statistics and various overhead functions
int main(int argc, char** argv)
{
    int num_threads = 1;
    int repetitions = 1;
    
    // parse input
    if(argc != 3)
    {
        printf("USAGE: hot <number of threads> <number of repetitions>\n");
        exit(1);
    }
    
    if(sscanf(argv[1], "%d", &num_threads) == EOF || num_threads < 1)
    {
        printf("Error: First argument must be a valid positive integer.\n");
        exit(1);
    }
    
    if(sscanf(argv[2], "%d", &repetitions) == EOF || repetitions < 1)
    {
        printf("Error: Second argument must be a valid positive integer.\n");
        exit(1);
    }
    
    // variables for statistics
    double total_time = 0;
    double fastest_time = FLT_MAX;
    
    int r;
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
        run_hotplate(num_threads, size, error, &iteration_count, &cell_count_gt_50_degrees);
    
        // get stop time
        double end_time = get_seconds();
        double time_interval = end_time - start_time;
        
        // report convergence and time
        printf("Iterations: %d\n", iteration_count);
        //printf("Cells with >= 50.0 degrees: %d\n", cell_count_gt_50_degrees);
        //printf("%d %f\n", num_threads, time_interval); // number_of_threads time_to_execute
        //fflush(stdout);
        
        total_time += time_interval;
        if(time_interval < fastest_time)
        {
            fastest_time = time_interval;
        }
    }
    
    // report average and fastest times
    printf("%d %f\n", num_threads, fastest_time);
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

void run_hotplate(int num_threads, int size, float error, int* iterations, int* cell_count_gt_50_degrees)
{
    
    //double start_time = get_seconds();
    
    (*iterations) = 0;
    (*cell_count_gt_50_degrees) = 0;
    
    int keep_going = 1;
    int rc = 0;
    void* status = &rc;
    
    // allocate space for threads
    ThreadInfo* threads_info = malloc(sizeof(ThreadInfo)*num_threads);
    pthread_t* threads = malloc(sizeof(pthread_t)*num_threads);
    // allocate space for plates
    float* current_plate = malloc(size*size*sizeof(float));
    float* next_plate = malloc(size*size*sizeof(float));
    int* test = malloc(size*size*sizeof(int));
    // allocate space for locks
    pthread_mutex_t* chunk_lock = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_t* count_lock = malloc(sizeof(pthread_mutex_t));
    
    
    // initialize the matrices
    initialize(size, (float(*) []) current_plate);
    initialize(size, (float(*) []) next_plate);
    initialize_test_cells(size, (int(*) []) test);
    // initialize the locks
    pthread_mutex_init(chunk_lock, NULL);
    pthread_mutex_init(count_lock, NULL);
    
    int count = 0;
    int row_index = 1;
    int col_index = 1;
    
    int chunk_index = 0;
    int chunk_size = 64;
    int chunk_count = size/chunk_size;
    if(size%chunk_size != 0)
    {
        chunk_count++;
    }
    
    // init thread info
    int i;
    for(i = 0; i < num_threads; i++)
    {
        threads_info[i].my_thread_num = i;
        threads_info[i].num_of_threads = num_threads;
        
        threads_info[i].size = size;
        threads_info[i].current_hotplate = &current_plate;
        threads_info[i].next_hotplate = &next_plate;
        threads_info[i].chunk_size = chunk_size;
        threads_info[i].chunk_index = &chunk_index;
        threads_info[i].chunk_count = chunk_count;
        threads_info[i].chunk_lock = chunk_lock;
        
        threads_info[i].error = error;
        threads_info[i].test_plate = &test;
        threads_info[i].keep_going = &keep_going;
        threads_info[i].row_index = &row_index;
        threads_info[i].col_index = &col_index;
        
        threads_info[i].iterations = iterations;
        threads_info[i].cell_count = cell_count_gt_50_degrees;
        threads_info[i].my_compute_time = 0.0f;
        threads_info[i].my_converged_time = 0.0f;
        threads_info[i].my_converged_count = 0;
        threads_info[i].my_wait_time = 0.0f;
        
        threads_info[i].count = &count;
        threads_info[i].count_lock = count_lock;
        
        threads_info[i].all_threads_info = threads_info;
    }
    
    // create threads
    for(i = 1; i < num_threads; i++)
    {
        rc = pthread_create(&threads[i], NULL, run_thread, (void*)&threads_info[i]);
        if(rc)
        {
            printf("ERROR: Return code from pthread_create() is %d.\n", rc);
            exit(-1);
        }
    }
    run_thread(&threads_info[0]);
    
    // join all threads
    for(i = 1; i < num_threads; i++)
    {
        rc = pthread_join(threads[i], &status);
        if(rc)
        {
            printf("ERROR: Return code from pthread_join() is %d.\n", rc);
            exit(-1);
        }
        //printf("Main: Joined thread %d with status of %ld.\n", i, (long) status);
    }
    
    // free allocated thread info
    free(threads_info);
    free(threads);
    // free allocated plates
    free(current_plate);
    free(next_plate);
    free(test);
    // destroy and free allocated locks
    pthread_mutex_destroy(chunk_lock);
    free(chunk_lock);
    pthread_mutex_destroy(count_lock);
    free(count_lock);
}

void* run_thread(void* thread_info)
{
    ThreadInfo* info = (ThreadInfo*) thread_info;
    //printf("Thread %d of %d running.\n", info->my_thread_num, info->num_of_threads);
    
    double start_time;
    int next_chunk_index = -1;
    
    while(*(info->keep_going))
    {
        // perform computations until nothing is left
        /*while((*(info->chunk_index)) != info->chunk_count)
        {
            pthread_mutex_lock(info->chunk_lock);
            if((*(info->chunk_index)) != info->chunk_count)
            {
                next_chunk_index = *(info->chunk_index);
                (*(info->chunk_index))++;
            }
            else
            {
                next_chunk_index = info->chunk_count;
            }
            pthread_mutex_unlock(info->chunk_lock);
            
            if(next_chunk_index != info->chunk_count)
            {
                double start_time = get_seconds();
                compute(next_chunk_index, info->chunk_size, info->size, (float(*)[])(*(info->current_hotplate)), (float(*)[])(*(info->next_hotplate)));
                info->my_compute_time += (get_seconds() - start_time);
            }
        }*/
        
        start_time = get_seconds();
        compute_static(info, info->size, (float(*)[])(*(info->current_hotplate)), (float(*)[])(*(info->next_hotplate)));
        info->my_compute_time += (get_seconds() - start_time);
        
        start_time = get_seconds();
        barrier(info);
        info->my_wait_time += (get_seconds() - start_time);
    }
    
    printf("Thread %d had %lf %lf %d %lf\n", info->my_thread_num, info->my_compute_time, info->my_converged_time, info->my_converged_count, info->my_wait_time);
    
    if((info->my_thread_num) != 0)
    {
        pthread_exit(NULL);
    }
    else
    {
        return 0;
    }
}

inline void swap(float** current, float** next)
{
    float* temp = *current;
    *current = *next;
    *next = temp;
}

void compute(int next_chunk_index, int chunk_size, int size, float current[size][size], float next[size][size])
{
    int row, col, row_end, col_end;
    row = chunk_size*next_chunk_index;
    row_end = row + chunk_size;
    col_end = size - 1;
    if(row < 1)
    {
        row = 1;
    }
    if(row_end > col_end)
    {
        row_end = col_end;
    }
    
    for(; row < row_end; row++)
    {
        float* top = current[row+1];
        float* curr = current[row];
        float* bottom = current[row-1];
        for(col = 1; col < col_end; col++)
        {
            next[row][col] = (bottom[col] + top[col] + curr[col - 1] + curr[col + 1] + 4.0f*curr[col])/8.0f;
        }
    }
}

int quick_convergence_test(ThreadInfo* info, int size, float plate[size][size], float error, int test[size][size])
{
    int row, col, converged, first_iteration;
    
    row = *(info->row_index);
    col = *(info->col_index);
    first_iteration = 1;
    converged = 1;
    
    for(; row < size - 1; row++)
    {
        if(!first_iteration)
        {
            col = 1;
        }
        else
        {
            first_iteration = 0;
        }
        
        for(; col < size - 1; col++)
        {
            if(!test[row][col])
            {
                //(info->my_converged_count)++;
                float average = (plate[row - 1][col] + plate[row + 1][col] + 
                                 plate[row][col - 1] + plate[row][col + 1])/4.0f;
                float difference = fabsf(plate[row][col] - average);
                    
                if(difference >= error)
                {
                    *(info->row_index) = row;
                    *(info->col_index) = col;
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

void barrier(ThreadInfo* info)
{
    info->my_wait = 1;
    
    // get the lock
    pthread_mutex_lock(info->count_lock);
    (*(info->count))++;
    
    // if I'm the last one in do the check
    if(*(info->count) == info->num_of_threads)
    {
        pthread_mutex_unlock(info->count_lock);
        
        double start_time = get_seconds();
        *(info->keep_going) = !quick_convergence_test(info, info->size, (float(*)[])(*(info->current_hotplate)), info->error, (int(*)[])(*(info->test_plate)));
        info->my_converged_time += (get_seconds() - start_time);
        
        *(info->count) = 0;
        
        if(*(info->keep_going))
        {
            swap(info->current_hotplate, info->next_hotplate);
            (*(info->iterations))++;
            //printf("Iter: %d\n", *(info->iterations));
        }
        else
        {
            *(info->row_index) = 1;
            *(info->col_index) = 1;
            *(info->keep_going) = !quick_convergence_test(info, info->size, (float(*)[])(*(info->current_hotplate)), info->error, (int(*)[])(*(info->test_plate)));
            if(*(info->keep_going))
            {
                swap(info->current_hotplate, info->next_hotplate);
                (*(info->iterations))++;
            }
        }
        
        int i;
        for(i = 0; i < info->num_of_threads; i++)
        {
            ((info->all_threads_info)[i]).my_wait = 0;
        }
    }
    else
    {
        pthread_mutex_unlock(info->count_lock);
        while(info->my_wait);// printf("Thread %d\n", info->my_thread_num);
    }
}

void compute_static(ThreadInfo* info, int size, float current[size][size], float next[size][size])
{
    int row, col, end, chunk_size;
    chunk_size = size/(info->num_of_threads);
    row = chunk_size*(info->my_thread_num);
    end = row + chunk_size;
    if(row < 1)
    {
        row = 1;
    }
    if(end > size - 1 || ((info->my_thread_num) == (info->num_of_threads - 1)))
    {
        end = size - 1;
    }
    
    //printf("Thread %d computing from [%d to %d)\n", info->my_thread_num, row, end);
    
    for(; row < end; row++)
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
