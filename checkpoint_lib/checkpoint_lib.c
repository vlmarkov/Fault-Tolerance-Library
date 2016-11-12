#include "checkpoint_lib.h"

double GLOBAL_START_TIME = 0.0;

double wtime_()
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

void make_snapshot(void *data, int count, MPI_Datatype datatype, int phase)
{
    MPI_File snapshot;
    MPI_Status status;

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
                   MPI_INFO_NULL,&snapshot );

    MPI_File_write(snapshot, data, count, datatype, &status);

    MPI_File_close(&snapshot);
}
