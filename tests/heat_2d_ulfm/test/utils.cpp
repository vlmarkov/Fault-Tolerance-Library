#include <iostream>

#include <mpi.h>
#include <stdlib.h>

void *xcalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);
    if (!p)
    {
        std::cerr << "No enough memory" << std::endl;
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
