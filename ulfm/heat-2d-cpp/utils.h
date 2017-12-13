#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdlib.h>

void *xcalloc(size_t nmemb, size_t size);

int get_block_size(int n, int rank, int nprocs);

int get_sum_of_prev_blocks(int n, int rank, int nprocs);

int check_overflow(int idx, int border);

#endif /* _UTILS_H_ */
