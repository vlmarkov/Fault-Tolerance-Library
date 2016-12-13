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
double cpl_start_time = 0.0;

void   **cpl_checkpoint_table;

int cpl_size    = 0;
int cpl_time    = 0;
int cpl_counter = 0;

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
    static int counter  = 0;
    char file_name[256] = { 0 };
    char file_path[256] = { 0 };
    /* 
     * 1_1 => [PHASE_OF_CALCULATION]_[COUNTER]
     *
     * PHASE_OF_CALCULATION - global state, need for 'goto'
     * COUNTER              - each PHASE_OF_CALCULATION could reach many times
     */

    sprintf(file_path,"%s/%d", SNAPSHOT_DIR_NAME, get_comm_rank__());
    mkdir(file_path, 0777);

    sprintf(file_name,"%s/%d/%d_%d", SNAPSHOT_DIR_NAME, get_comm_rank__(), phase, counter++);

    MPI_File_open( MPI_COMM_WORLD, file_name, 
                   MPI_MODE_CREATE|MPI_MODE_WRONLY, 
                   MPI_INFO_NULL, snapshot );
}

void close_snapshot_file_(MPI_File *snapshot)
{
    MPI_File_close(snapshot);
}

void write_to_snapshot_(MPI_File file, void *data, int n, MPI_Datatype type)
{
    MPI_Status status;
    MPI_File_write(file, data, n, type, &status);
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
        //printf("%s\n", line);
        char tmp[10] = { 0 };
        sprintf(tmp, "%c", line[0]);
        if (atoi(tmp) == myrank) {
            //printf("*%s\n", line);
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

    printf("Rankd %d, file %s, phase %d\n", myrank, last_checkpoint, last_checkpoint[0] - '0');

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