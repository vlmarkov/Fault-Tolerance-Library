#include "utils.h"

#include <mpi.h>
#include <stdio.h>

void *xcalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);
    if (!p)
    {
        fprintf(stderr, "No enough memory\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    return p;
}

int get_block_size(int n, int rank, int nprocs)
{
    int s = n / nprocs;

    if (n % nprocs > rank)
    {
        s++;
    }

    return s;
}

int get_sum_of_prev_blocks(int n, int rank, int nprocs)
{
    int rem = n % nprocs;

    return n / nprocs * rank + ((rank >= rem) ? rem : rank);
}
