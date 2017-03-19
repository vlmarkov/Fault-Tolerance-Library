/*****************************************************************************/
/* C - check                                                                 */
/* P - point                                                                 */
/* L - library                                                               */
/*****************************************************************************/

#include "checkpoint_lib.h"

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>

#include <sys/stat.h>

/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
enum {
    CPL_CHECKPOINT_MODE = 1,
    CPL_RECOVERY_MODE   = 2
};


void   **cpl_checkpoint_table;

double cpl_start_time    = 0.0;
double cpl_start_time_local = 0.0;

double cpl_save_time       = 0.0;
double cpl_save_time_local = 0.0;

int cpl_snapshot_counter = 0;

int cpl_size    = 0;
int cpl_time    = 0;
int cpl_counter = 0;

int cpl_run_options = 0;

static struct itimerval cpl_timer;


/*****************************************************************************/
/* Timer                                                                     */
/*****************************************************************************/
inline void timer_init_()
{
    cpl_timer.it_interval.tv_sec  = cpl_time; // interval 
    cpl_timer.it_interval.tv_usec = 0;
    cpl_timer.it_value.tv_sec     = cpl_time; // time until next expiration
    cpl_timer.it_value.tv_usec    = 0;

    setitimer(ITIMER_REAL, &cpl_timer, NULL);
}

inline void timer_stop_()
{
    cpl_timer.it_interval.tv_sec = 0;
    cpl_timer.it_value.tv_sec    = 0;
}


/*****************************************************************************/
/* Measure time                                                              */
/*****************************************************************************/
inline double wtime_()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}


/*****************************************************************************/
/* Aditional internal mpi functions                                          */
/*****************************************************************************/
static int get_comm_rank__()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

static int get_comm_size__()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}


/*****************************************************************************/
/* Work with files                                                           */
/*****************************************************************************/
void open_snapshot_file_(MPI_File *snapshot, int phase)
{
    MPI_Status status;

    static int counter  = 0;
    char file_name[256] = { 0 };
    char file_path[256] = { 0 };

    /*
     * 1_1 => [PHASE_OF_CALCULATION]_[COUNTER]
     *
     * PHASE_OF_CALCULATION - global state, need for 'goto'
     * COUNTER              - each PHASE_OF_CALCULATION could reach many times
     */

    sprintf(file_path,"/tmp/%s/%d", SNAPSHOT_DIR_NAME, get_comm_rank__());
    mkdir(file_path, 0777);

    sprintf(file_name,"%s/%d/%d_%d", SNAPSHOT_DIR_NAME, get_comm_rank__(), phase, counter++);

    MPI_File_open( MPI_COMM_WORLD, file_name, 
                   MPI_MODE_CREATE|MPI_MODE_WRONLY, 
                   MPI_INFO_NULL, snapshot );

    cpl_save_time_local = wtime_();
}

void close_snapshot_file_(MPI_File *snapshot)
{
    MPI_Status status;

    cpl_snapshot_counter += 1;

    cpl_save_time_local = wtime_() - cpl_save_time_local;
    cpl_save_time += cpl_save_time_local;

    double elapsed_time = wtime_() - cpl_start_time_local;

    MPI_File_write(*snapshot, &elapsed_time, 1, MPI_DOUBLE, &status);
    MPI_File_write(*snapshot, &cpl_save_time, 1, MPI_DOUBLE, &status);
    MPI_File_write(*snapshot, &cpl_snapshot_counter, 1, MPI_INT, &status);
    MPI_File_write(*snapshot, &INTEGRITY_SNAPSHOT, strlen(INTEGRITY_SNAPSHOT), MPI_CHAR, &status);

    MPI_File_close(snapshot);
}

void write_to_snapshot_(MPI_File file, void *data, int n, MPI_Datatype type)
{
    MPI_Status status;
    MPI_File_write(file, data, n, type, &status);
}

void open_snapshot_file_tmp_(FILE **snapshot, int phase)
{
    static int counter  = 0;
    char file_name[256] = { 0 };
    char file_path[256] = { 0 };

    /*
     * 1_1_1 => [RANK]_[PHASE_OF_CALCULATION]_[COUNTER]
     *
     * PHASE_OF_CALCULATION - global state, need for 'goto'
     * COUNTER              - each PHASE_OF_CALCULATION could reach many times
     */
    sprintf(file_name,"/tmp/snapshot_%d_%d_%d", get_comm_rank__(), phase, counter++);

    *snapshot = fopen(file_name, "w");
    if (!*snapshot) {
        fprintf(stderr, "Can't create snapshot %d\n", get_comm_rank__());
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    cpl_save_time_local = wtime_();
}

void close_snapshot_file_tmp_(FILE *snapshot)
{
    cpl_snapshot_counter += 1;

    cpl_save_time_local = wtime_() - cpl_save_time_local;
    cpl_save_time += cpl_save_time_local;

    double elapsed_time = wtime_() - cpl_start_time_local;

    fprintf(snapshot, "%f\n", elapsed_time);
    fprintf(snapshot, "%f\n", cpl_save_time);
    fprintf(snapshot, "%d\n", cpl_snapshot_counter);
    fprintf(snapshot, "%s\n", INTEGRITY_SNAPSHOT);

    fclose(snapshot);
}

void write_to_snapshot_tmp_(FILE *snapshot, void *data, int n, MPI_Datatype type)
{
    switch (type)
    {
        case MPI_DOUBLE:
            fwrite(data, sizeof(double), n, snapshot);
            break;
        case MPI_INT:
            fwrite(data, sizeof(int), n, snapshot);
            break;
        case MPI_CHAR:
            fwrite(data, sizeof(char), n, snapshot);
            break;
        default:
            fprintf(stderr, "Can't write to file - unsupported type %d\n", type);
            break;
    }
}

FILE *cpl_open_file(char *file_name, char *mode)
{
    FILE *file = fopen(file_name, mode);
    if (!file) {
        fprintf(stderr, "Can't read %s\n", file_name);
        exit(1);
    }

    long offset = sizeof(cpl_start_time) \
                  + sizeof(cpl_save_time) \
                  + sizeof(cpl_snapshot_counter) \
                  + strlen(INTEGRITY_SNAPSHOT); 

    // Move file to the end of file
    fseek(file, -offset, SEEK_END);

    // copy the file into the buffer:
    fread(&cpl_start_time, sizeof(double), 1, file);
    fread(&cpl_save_time, sizeof(double), 1, file);
    fread(&cpl_snapshot_counter, sizeof(int), 1, file);

    // Move file ptr to the begin of file
    rewind(file);

    return file;
}

/*****************************************************************************/
/* Get last snapshot file                                                    */
/*****************************************************************************/
int get_last_snapshot_(char *last_checkpoint)
{
    int myrank = get_comm_rank__();

    char tmp_checkpoint[256] = { 0 };

    char * line = NULL;
    size_t len = 0;
    int i, j;

    FILE *file = fopen(INTEGRITY_SNAPSHOT_FILE, "r");
    if (!file) {
        fprintf(stderr, "Can't read from %s\n", INTEGRITY_SNAPSHOT_FILE);
        exit(1);
    }

    while ((getline(&line, &len, file)) != -1) {
        char tmp[10] = { 0 };
        sprintf(tmp, "%c", line[0]);
        if (atoi(tmp) == myrank) {
            strcpy(last_checkpoint, line);
        }
    }
    fclose(file);

    for (i = 0; i < strlen(last_checkpoint); i++) {
        if (last_checkpoint[i] == '=') {
            break;
        }
    }

    i += 1;

    for (j = 0; i < strlen(last_checkpoint); j++, i++) {
        tmp_checkpoint[j] = last_checkpoint[i];
    }

    strcpy(last_checkpoint, tmp_checkpoint);

    last_checkpoint[strlen(last_checkpoint) - 1] = '\0';

    printf("Rank %d, file %s, phase %d\n", myrank, last_checkpoint, last_checkpoint[0] - '0');

    return last_checkpoint[0] - '0';
}


/*****************************************************************************/
/* Get checkpoint index by checkpoint name                                   */
/*****************************************************************************/
int get_checkpoint_idx_by_name_(void **table, int size, void *name)
{
    int i;
    for (i = 0; i < size; i++) {
        if (table[i] == name) {
            break;
        }
    }
    return i;
}


/*****************************************************************************/
/* Init                                                                      */
/*****************************************************************************/
void **init_table_(int size)
{
    void **jump_table = (void**) malloc (sizeof(void*) * size);
    if (!jump_table) {
        fprintf(stderr, "[ERROR] can't allocate memory for CPL_GLOBAL_JUMP_TABLE\n");
        exit(1);
    }

    mkdir(SNAPSHOT_DIR_NAME, 0777);

    return jump_table;
}

void cpl_init(int size, double time, int argc, char *argv[])
{
    cpl_size             = time; 
    cpl_counter          = 0;
    cpl_size             = size;
    cpl_start_time_local = wtime_();
    cpl_checkpoint_table = init_table_(cpl_size);

    if (get_comm_rank__() == 0) {
        printf("\n[CPL_LIBRARY] note, use: %s [args] [recovery]\n", argv[0]);
    }

    while (argc-->0) {
        if (strcmp(argv[argc], "recovery") == 0) {
            if (get_comm_rank__() == 0) {
                printf("[CPL_LIBRARY] running options 'recovery'\n");
            }
            cpl_run_options = CPL_RECOVERY_MODE;
            break;
        }
    }
}

void cpl_finalize()
{
    int i;

    // Zero process print statistic
    if (get_comm_rank__() == 0) {
        printf("\n");
        // Print up border
        for (i = 0; i < 80; i++) {
            printf("*");
        }
        printf("\n");
        cpl_start_time +=  wtime_() - cpl_start_time_local;
        printf("[CPL_LIBRARY] common elapsed time   = %f sec\n", cpl_start_time);

        double time_one_save = cpl_save_time / (double)cpl_snapshot_counter;
        cpl_save_time = time_one_save * (double)cpl_snapshot_counter;

        printf("[CPL_LIBRARY] elapsed time one save = %f sec\n", time_one_save);
        printf("[CPL_LIBRARY] elapsed time all save = %f sec\n", cpl_save_time);
        printf("[CPL_LIBRARY] amount of snapshots   = %d\n", cpl_snapshot_counter);

        // Print down border
        for (i = 0; i < 80; i++) {
            printf("*");
        }
        printf("\n");
    }
}

int IS_CPL_RECOVERY_MODE()
{
    if (cpl_run_options == CPL_RECOVERY_MODE)
        return 1;
    else
        return 0;
}
