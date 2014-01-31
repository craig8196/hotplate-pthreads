/*
 * Each 2D array is organized such that the first index is the row, the second is the column.
 * 
 * Version 4: Busy wait.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <float.h>

// structs
typedef struct linear_barrier_t
{
    int count;
    pthread_mutex_t* count_lock;
    pthread_cond_t* count_cond;
} LinearBarrier;

typedef struct thread_info_t
{
    int my_thread_num;
    int num_of_threads;
    int size;
    float error;
    float** current_hotplate;
    float** next_hotplate;
    int** test_plate;
    int* keep_going;
    int my_keep_going;
    int* iterations;
    int* cell_count;
    double my_compute_time;
    double my_converged_time;
    int my_converged_count;
    double my_wait_time;
    struct thread_info_t* all_threads_info;
    LinearBarrier* barrier;
    volatile int my_wait;
} ThreadInfo;



// function signatures
double get_seconds();
void initialize(int size, float plate[size][size]);
void initialize_test_cells(int size, int test[size][size]);
void set_static_cells(int size, float plate[size][size]);
void swap(float** current, float** next);
void compute(ThreadInfo* info, int size, float current[size][size], float next[size][size]);
void has_converged(ThreadInfo* info, int size, float plate[size][size], float error, int test[size][size]);
void* run_thread(void* thread_info);
void run_hotplate(int num_threads, int size_of_plate, float error, int* iterations, int* cell_count_gt_50_degrees);
void barrier(ThreadInfo* info);
LinearBarrier* new_linear_barrier();
void destroy_linear_barrier(LinearBarrier** barrier);


// performs statistics and various overhead functions
int main(int argc, char** argv)
{
    int num_threads = 1;
    int repetitions = 1;
    
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
        printf("Cells with >= 50.0 degrees: %d\n", cell_count_gt_50_degrees);
        printf("%d %f\n", num_threads, time_interval); // number_of_threads time_to_execute
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
    
    double start_time = get_seconds();
    
    (*iterations) = 0;
    (*cell_count_gt_50_degrees) = 0;
    
    int keep_going = 1;
    int rc = 0;
    void* status = &rc;
    
    // allocate space for threads
    ThreadInfo* threads_info = malloc(sizeof(ThreadInfo)*num_threads);
    pthread_t* threads = malloc(sizeof(pthread_t)*num_threads);
    // allocate space for
    float* current_plate = malloc(size*size*sizeof(float));
    float* next_plate = malloc(size*size*sizeof(float));
    int* test = malloc(size*size*sizeof(int));
    LinearBarrier* barrier = new_linear_barrier(barrier);
    
    
    // initialize the matrices
    initialize(size, (float(*) []) current_plate);
    initialize(size, (float(*) []) next_plate);
    initialize_test_cells(size, (int(*) []) test);
    
    
    // init thread info
    int i;
    for(i = 0; i < num_threads; i++)
    {
        threads_info[i].my_thread_num = i;
        threads_info[i].num_of_threads = num_threads;
        threads_info[i].size = size;
        threads_info[i].error = error;
        threads_info[i].current_hotplate = &current_plate;
        threads_info[i].next_hotplate = &next_plate;
        threads_info[i].test_plate = &test;
        threads_info[i].keep_going = &keep_going;
        threads_info[i].my_keep_going = 1;
        threads_info[i].iterations = iterations;
        threads_info[i].cell_count = cell_count_gt_50_degrees;
        threads_info[i].my_compute_time = 0.0f;
        threads_info[i].my_converged_time = 0.0f;
        threads_info[i].my_converged_count = 0;
        threads_info[i].my_wait_time = 0.0f;
        threads_info[i].all_threads_info = threads_info;
        threads_info[i].barrier = barrier;
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
    
    double setup_time = get_seconds() - start_time;
    
    run_thread(&threads_info[0]);
    
    start_time = get_seconds();
    
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
    
    //printf("Iter: %d\n", *iterations);
    
    // free allocated memory
    free(threads_info);
    free(threads);
    free(current_plate);
    free(next_plate);
    free(test);
    destroy_linear_barrier(&barrier);
    
    double teardown_time = get_seconds() - start_time;
    
    printf("Setup and teardown taken about %lf seconds.\n", (teardown_time + setup_time));
}

void* run_thread(void* thread_info)
{
    ThreadInfo* info = (ThreadInfo*) thread_info;
    //printf("Thread %d of %d running.\n", info->my_thread_num, info->num_of_threads);
    
    
    while(*(info->keep_going))
    {
        double start_time = get_seconds();
        compute(info, info->size, (float(*)[])(*(info->current_hotplate)), (float(*)[])(*(info->next_hotplate)));
        info->my_compute_time += (get_seconds() - start_time);
        start_time = get_seconds();
        has_converged(info, info->size, (float(*)[])(*(info->current_hotplate)), info->error, (int(*)[])(*(info->test_plate)));
        info->my_converged_time += (get_seconds() - start_time);
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

void compute(ThreadInfo* info, int size, float current[size][size], float next[size][size])
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

void has_converged(ThreadInfo* info, int size, float plate[size][size], float error, int test[size][size])
{
    int row, col, end, chunk_size;
    info->my_keep_going = 0;
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
    
    for(; row < end; row++)
    {
        for(col = 1; col < size - 1; col++)
        {
            if(!test[row][col])
            {
                (info->my_converged_count)++;
                float average = (plate[row - 1][col] + plate[row + 1][col] + 
                                 plate[row][col - 1] + plate[row][col + 1])/4.0f;
                           
                
                float difference = fabsf(plate[row][col] - average);
                
                /*printf("Avg: %f", average);
                printf("Val: %f", plate[row][col]);
                printf("Dif: %f", difference);*/
                    
                if(difference >= error)
                {
                    info->my_keep_going = 1;
                    break;
                }
            }
        }
        if(info->my_keep_going)
        {
            break;
        }
    }
}

void barrier(ThreadInfo* info)
{
    LinearBarrier* barrier = info->barrier;
    
    info->my_wait = 1;
    pthread_mutex_lock(barrier->count_lock);
    (barrier->count)++;
    
    if(barrier->count == info->num_of_threads)
    {
        pthread_mutex_unlock(barrier->count_lock);
        
        barrier->count = 0;
        *(info->keep_going) = 0;
        
        int i;
        for(i = 0; i < info->num_of_threads; i++)
        {
            *(info->keep_going) = *(info->keep_going) || ((info->all_threads_info)[i]).my_keep_going;
        }
        
        if(*(info->keep_going))
        {
            swap(info->current_hotplate, info->next_hotplate);
            (*(info->iterations))++;
        }
        
        for(i = 0; i < info->num_of_threads; i++)
        {
            ((info->all_threads_info)[i]).my_wait = 0;
        }
    }
    else
    {
        pthread_mutex_unlock(barrier->count_lock);
        while(info->my_wait){}// printf("Thread %d\n", info->my_thread_num);
    }
}


LinearBarrier* new_linear_barrier()
{
    LinearBarrier* barrier = malloc(sizeof(LinearBarrier));
    
    barrier->count_lock = malloc(sizeof(pthread_mutex_t));
    barrier->count_cond = malloc(sizeof(pthread_cond_t));
    barrier->count = 0;
    pthread_mutex_init(barrier->count_lock, NULL);
    pthread_cond_init(barrier->count_cond, NULL);
    
    return barrier;
}
void destroy_linear_barrier(LinearBarrier** barrier)
{
    pthread_mutex_destroy((*barrier)->count_lock);
    pthread_cond_destroy((*barrier)->count_cond);
    
    free((*barrier)->count_lock);
    free((*barrier)->count_cond);
    free(*barrier);
    (*barrier) = NULL;
}
