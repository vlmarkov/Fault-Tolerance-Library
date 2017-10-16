#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdlib.h>

#define EPS 0.001
#define PI 3.14159265358979323846

#define NELEMS(x) (sizeof((x)) / sizeof((x)[0]))
#define IND(i, j) ((i) * (nx + 2) + (j))

void *xcalloc(size_t nmemb, size_t size);

int get_block_size(int n, int rank, int nprocs);

int get_sum_of_prev_blocks(int n, int rank, int nprocs);

#endif /* _UTILS_H_ */
