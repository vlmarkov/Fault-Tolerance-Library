#include "Utils.h"

#ifdef MPI_SUPPORT
#include <mpi.h>
#endif /* MPI_SUPPORT */

#include <iostream>

void *xCalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);

    if (!p)
    {
        std::cerr << "<" << __FUNCTION__ << ">"
                  << " No enough memory"
                  << std::endl;

#ifdef MPI_SUPPORT
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

    }

    return p;
}

int getBlockSize(int n, int rank, int nprocs)
{
    int s = n / nprocs;

    if (n % nprocs > rank)
    {
        s++;
    }

    return s;
}

int getSumOfPrevBlocks(int n, int rank, int nprocs)
{
    int rem = n % nprocs;
    return n / nprocs * rank + ((rank >= rem) ? rem : rank);
}

int checkOverflow(int idx, int border)
{
    if (idx > border - 1)
    {
        idx = idx - border;
    }

    return idx;
}
