/*
 * Each 2D array is organized such that the first index is the row, the second is the column.
 * 
 * Version 5: Log barrier.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <float.h>

// structs
typedef struct barrier_node {
        pthread_mutex_t count_lock;
        pthread_cond_t ok_to_proceed_up;
        pthread_cond_t ok_to_proceed_down;
        int count;
} BarrierNode;

typedef struct log_barrier {
    BarrierNode* barrier_nodes;
    int number_in_barrier;
    pthread_mutex_t logbarrier_count_lock;
} LogBarrier;

typedef struct thread_info_t
{
    int my_thread_num;
    int num_of_threads;
    int size;
    float error;
    float** current_hotplate;
    float** next_hotplate;
    int** test_plate;
    volatile int* keep_going;
    volatile int my_keep_going;
    volatile int* iterations;
    int* cell_count;
    double my_compute_time;
    double my_converged_time;
    int my_converged_count;
    double my_wait_time;
    struct thread_info_t* all_threads_info;
    LogBarrier* barrier;
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
LogBarrier* new_log_barrier(int num_threads);
void destroy_log_barrier(LogBarrier** barrier);
void mylib_init_barrier_nodes(int num_threads, BarrierNode* b);



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
    
    volatile int keep_going = 1;
    int rc = 0;
    void* status = &rc;
    
    // allocate space for threads
    ThreadInfo* threads_info = malloc(sizeof(ThreadInfo)*num_threads);
    pthread_t* threads = malloc(sizeof(pthread_t)*num_threads);
    // allocate space for
    float* current_plate = malloc(size*size*sizeof(float));
    float* next_plate = malloc(size*size*sizeof(float));
    int* test = malloc(size*size*sizeof(int));
    LogBarrier* barrier = new_log_barrier(num_threads);
    
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
    destroy_log_barrier(&barrier);
    
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


LogBarrier* new_log_barrier(int num_threads)
{
    LogBarrier* barrier = malloc(sizeof(LogBarrier));
    
    barrier->barrier_nodes = malloc(sizeof(BarrierNode)*num_threads);
    
    mylib_init_barrier_nodes(num_threads, barrier->barrier_nodes);
    pthread_mutex_init(&(barrier->logbarrier_count_lock), NULL);
    barrier->number_in_barrier = 0;
    
    return barrier;
}
void destroy_log_barrier(LogBarrier** barrier)
{
    pthread_mutex_destroy(&((*barrier)->logbarrier_count_lock));
    
    free((*barrier)->barrier_nodes);
    free(*barrier);
    (*barrier) = NULL;
}

void mylib_init_barrier_nodes(int num_threads, BarrierNode* b)
{
        int i;
        for (i = 0; i < num_threads; i++) {
                b[i].count = 0;
                pthread_mutex_init(&(b[i].count_lock), NULL);
                pthread_cond_init(&(b[i].ok_to_proceed_up), NULL);
                pthread_cond_init(&(b[i].ok_to_proceed_down), NULL);
        }
}

void barrier (ThreadInfo* info) //mylob_logbarrier_t b, int num_threads, int thread_id)
{
    LogBarrier* barrier = info->barrier;
    BarrierNode* b = barrier->barrier_nodes;
    
    int num, base, index;
    num = 2;
    base = 0;

    if (info->num_of_threads == 1)
    {
        *(info->keep_going) = 0;
        
        int i;
        for(i = 0; i < info->num_of_threads; i++)
        {
            *(info->keep_going) += (info->all_threads_info)[i].my_keep_going;
        }
        
        if(*(info->keep_going))
        {
            swap(info->current_hotplate, info->next_hotplate);
            (*(info->iterations))++;
            //printf("Iter: %d", *(info->iterations));
        }
        return;
    }
    
    pthread_mutex_lock(&(barrier->logbarrier_count_lock));
    (barrier->number_in_barrier)++;
    
    if (barrier->number_in_barrier == info->num_of_threads)
    {
        *(info->keep_going) = 0;
        (barrier->number_in_barrier) = 0;
        int i;
        for(i = 0; i < info->num_of_threads; i++)
        {
            *(info->keep_going) += (info->all_threads_info)[i].my_keep_going;
        }
        
        if(*(info->keep_going))
        {
            swap(info->current_hotplate, info->next_hotplate);
            (*(info->iterations))++;
            //printf("Iter: %d", *(info->iterations));
        }
    }
    pthread_mutex_unlock(&(barrier->logbarrier_count_lock));
    
    //printf("Thread %d\n", info->my_thread_num);
    
    do
    {
        index = base + (info->my_thread_num) /num;
        
        if ((info->my_thread_num) % num == 0)
        {
            pthread_mutex_lock(&(b[index].count_lock));
            b[index].count ++;
            while (b[index].count < 2)
                pthread_cond_wait(&(b[index].ok_to_proceed_up),
                                  &(b[index].count_lock));
                pthread_mutex_unlock(&(b[index].count_lock));
        }
        else
        {
            pthread_mutex_lock(&(b[index].count_lock));
            b[index].count ++;
            if (b[index].count == 2)
                pthread_cond_signal(&(b[index].ok_to_proceed_up));

            while (pthread_cond_wait(&(b[index].ok_to_proceed_down),
                                     &(b[index].count_lock)) != 0);
            pthread_mutex_unlock(&(b[index].count_lock));
            break;
        }
        base = base + info->num_of_threads/num;
        num =num * 2;
    } while (num <= info->num_of_threads);

    num = num / 2;

    for (; num > 1; num = num / 2)
    {
        base = base - info->num_of_threads/num;
        index = base + info->my_thread_num / num;
        pthread_mutex_lock(&(b[index].count_lock));
        b[index].count = 0;
        pthread_cond_signal(&(b[index].ok_to_proceed_down));
        pthread_mutex_unlock(&(b[index].count_lock));
    }
}
