#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdlib.h>

void *xCalloc(size_t nmemb, size_t size);

int getBlockSize(int n, int rank, int nprocs);

int getSumOfPrevBlocks(int n, int rank, int nprocs);

int checkOverflow(int idx, int border);

#endif /* _UTILS_H_ */
