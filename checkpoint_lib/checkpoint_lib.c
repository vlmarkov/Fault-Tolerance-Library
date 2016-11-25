#include "checkpoint_lib.h"

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>


double GLOBAL_START_TIME = 0.0;
struct itimerval nval, oval;

inline void timer_init_()
{
    nval.it_interval.tv_sec  = TIME; // interval 
    nval.it_interval.tv_usec = 0;
    nval.it_value.tv_sec     = TIME; // time until next expiration
    nval.it_value.tv_usec    = 0;

    setitimer(ITIMER_REAL, &nval, &oval);
}

inline void timer_stop_()
{
    nval.it_interval.tv_sec  = 0; // interval 
    nval.it_interval.tv_usec = 0;
    nval.it_value.tv_sec     = 0; // time until next expiration
    nval.it_value.tv_usec    = 0;
}

inline double wtime_()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

static int get_comm_rank_()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

static int get_comm_size_()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void open_checkpoint_file(MPI_File *snapshot, int phase)
{
    char file_name[256] = { 0 };

    /* 
     * 1_1_1.3456 => [PHASE_OF_CALCULATION]_[RANK]_[CHECKPOINT_TIME]
     *
     * PHASE_OF_CALCULATION - global state, need for 'goto'
     * RANK                 - process rank
     * CHECKPOINT_TIME      - each 'PHASE_OF_CALCULATION' could reach many times
     */

    sprintf(file_name,"%d_%d_%f", phase, get_comm_rank_(), wtime_() - GLOBAL_START_TIME);

    MPI_File_open( MPI_COMM_WORLD, file_name, 
                   MPI_MODE_CREATE|MPI_MODE_WRONLY, 
                   MPI_INFO_NULL, snapshot );
}

void close_checkpoint_file(MPI_File *snapshot)
{
    MPI_File_close(snapshot);
}

void make_snapshot(MPI_File file, void *data, int n, MPI_Datatype type)
{
    MPI_Status status;
    MPI_File_write(file, data, n, type, &status);
}

/*****************************************************************************/
/* Get last checkpoint file                                                  */
/*****************************************************************************/
int get_lastcheckpoint_rank_(char *file)
{
    int j, i, checkpoint_rank;

    char tmp_rank[10] = { 0 };

    for (j = 0, i = 2; i < strlen(file); i++, j++) {
        if (file[i] != '_') {
           tmp_rank[j] = file[i];
        } else {
            break;
        }
    }

    sscanf(tmp_rank, "%d", &checkpoint_rank);
    return checkpoint_rank;
}

void get_lastcheckpoint_time_(char *a, char *b)
{
    int j, i;

    double x, y;

    char time_last[10] = { 0 };
    char time_cur[10] = { 0 };

    for (j = 0, i = 4; i < strlen(b); i++, j++) {
        time_last[j] = a[i];
        time_cur[j]  = b[i];
    }
    
    sscanf(time_last, "%lf", &x);
    sscanf(time_cur, "%lf", &y);

    // Figure out which checkpoint is 'older'
    if (x < y) {
        strcpy(a, b);
    }
}

int get_lastcheckpoint(char *last_checkpoint)
{
    int myrank = get_comm_rank_();

    DIR           *dir;
    struct dirent *file;

    dir = opendir("./"); // open current directory
    if (dir) {
        while (file = readdir(dir)) {
            
            // Skip '.' '..' directory
            if (strlen(file->d_name) < 3) {
                continue;
            }

            if (file->d_name[1] != '_') {
                continue;
            }

            // Work only with myrank checkpoint
            if (myrank == get_lastcheckpoint_rank_(file->d_name)) {
                //printf("[DEBUG] Found for rank %d - %s\n", myrank, file->d_name);

                // Work only with greater checkpoint phase
                if (last_checkpoint[0] <= file->d_name[0]) {
                    get_lastcheckpoint_time_(last_checkpoint, file->d_name);
                }
                //printf("[DEBUG] %s\n", last_checkpoint);
            }
        }
        closedir(dir);
    } else {
        fprintf(stderr, "can't open current directory\n");
    }

    return last_checkpoint[0] - '0';
}



int get_checkpoint_idx_by_name(void **table, int size, void *name)
{
    int i;
    for (i = 0; i < size; i++) {
        if (table[i] == &name) {
            break;
        }
    }
}
