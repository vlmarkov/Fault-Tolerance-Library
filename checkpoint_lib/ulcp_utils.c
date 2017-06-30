#include "ulcp_utils.h"

#include <mpi.h>
#include <mpi-ext.h>

#include <sys/time.h>

double ulcp_wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

int ulcp_get_comm_rank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int ulcp_get_comm_size()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}
